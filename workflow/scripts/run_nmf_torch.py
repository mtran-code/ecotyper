import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch


EPSILON = 1e-8


def resolve_device(requested):
    requested = requested.lower()
    if requested == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")

    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("Requested CUDA NMF backend, but CUDA is not available")
    if requested == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("Requested MPS NMF backend, but MPS is not available")
    return torch.device(requested)


def read_posneg_matrix(input_path):
    raw = pd.read_csv(input_path, sep="\t", index_col=0)
    values = raw.to_numpy(dtype=np.float32, copy=True)
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    posneg = np.vstack([np.maximum(values, 0), np.maximum(-values, 0)])
    row_names = [f"{name}__pos" for name in raw.index] + [
        f"{name}__neg" for name in raw.index
    ]
    return posneg, row_names, list(raw.columns)


def kl_loss(v, w, h):
    reconstruction = torch.clamp(w @ h, min=EPSILON)
    return torch.sum(
        v * torch.log(torch.clamp(v, min=EPSILON) / reconstruction) - v + reconstruction
    )


def factorize(v_np, rank, seed, device, max_iter, tolerance):
    torch.manual_seed(seed)
    v = torch.tensor(v_np, dtype=torch.float32, device=device)
    rows, columns = v.shape
    w = torch.rand((rows, rank), dtype=torch.float32, device=device) + EPSILON
    h = torch.rand((rank, columns), dtype=torch.float32, device=device) + EPSILON

    previous = None
    final_loss = None
    iterations = 0
    for iteration in range(1, max_iter + 1):
        reconstruction = torch.clamp(w @ h, min=EPSILON)
        h = h * (
            (w.T @ (v / reconstruction))
            / torch.clamp(w.sum(dim=0)[:, None], min=EPSILON)
        )

        reconstruction = torch.clamp(w @ h, min=EPSILON)
        w = w * (
            ((v / reconstruction) @ h.T)
            / torch.clamp(h.sum(dim=1)[None, :], min=EPSILON)
        )

        if iteration == 1 or iteration % 10 == 0 or iteration == max_iter:
            current = float(kl_loss(v, w, h).detach().cpu())
            final_loss = current
            if previous is not None:
                relative_change = abs(previous - current) / max(abs(previous), EPSILON)
                if relative_change < tolerance:
                    iterations = iteration
                    break
            previous = current
        iterations = iteration

    return w.detach().cpu().numpy(), h.detach().cpu().numpy(), iterations, final_loss


def write_matrix(path, values, rows, columns):
    frame = pd.DataFrame(values, index=rows, columns=columns)
    frame.to_csv(path, sep="\t")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--rank", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--max-iter", default=500, type=int)
    parser.add_argument("--tolerance", default=1e-4, type=float)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    device = resolve_device(args.device)
    matrix, row_names, column_names = read_posneg_matrix(args.input)

    start = time.time()
    w, h, iterations, final_loss = factorize(
        matrix,
        args.rank,
        args.seed,
        device,
        args.max_iter,
        args.tolerance,
    )
    elapsed = time.time() - start

    component_names = [f"factor_{idx + 1}" for idx in range(args.rank)]
    write_matrix(output_dir / "torch_W.txt", w, row_names, component_names)
    write_matrix(output_dir / "torch_H.txt", h, component_names, column_names)
    with open(output_dir / "torch_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "backend": "torch",
                "device": str(device),
                "iterations": iterations,
                "final_loss": final_loss,
                "elapsed_seconds": elapsed,
                "rank": args.rank,
                "seed": args.seed,
            },
            handle,
            indent=2,
        )
        handle.write("\n")


if __name__ == "__main__":
    main()
