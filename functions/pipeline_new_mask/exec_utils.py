from __future__ import annotations

import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable


def run_command(cmd: str, check: bool = True, timeout: int | None = None, verbose: bool = True):
    """Run a shell command and return CompletedProcess."""
    if verbose:
        print(f"\n>>> {cmd}\n")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if result.stdout.strip() and verbose:
        print(result.stdout.strip())
    if result.stderr.strip() and verbose:
        print("STDERR:", result.stderr.strip())
    if check and result.returncode != 0:
        print(f"⚠️ Command failed with exit code {result.returncode}")
    return result


def check_file(path, label: str = "") -> bool:
    """Print existence/size summary for a file path and return status."""
    p = Path(path)
    if p.exists():
        size_mb = p.stat().st_size / (1024 * 1024)
        print(f"  ✅ {label or p.name}: {size_mb:.1f} MB")
        return True
    print(f"  ❌ MISSING: {label or p.name}")
    return False


def launch_gui(cmd: str) -> None:
    """Launch a GUI command in background without blocking."""
    print(f">>> {cmd}")
    subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _run_subject_cmd(args):
    """Worker for parallel execution with per-subject logs."""
    sid, cmd, timeout, log_dir, label = args
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        from datetime import datetime

        log_file = Path(log_dir) / f"{label}_{sid}.log"
        with open(log_file, "w", encoding="utf-8") as handle:
            handle.write(f"{'=' * 70}\n")
            handle.write(f"Subject: {sid}\n")
            handle.write(f"Label:   {label}\n")
            handle.write(f"Time:    {datetime.now().isoformat()}\n")
            handle.write(f"RC:      {result.returncode}\n")
            handle.write(f"Command: {cmd[:500]}\n")
            handle.write(f"{'=' * 70}\n\n")
            handle.write("=== STDOUT ===\n")
            handle.write(result.stdout)
            handle.write("\n\n=== STDERR ===\n")
            handle.write(result.stderr)

        return (
            sid,
            result.returncode,
            result.stdout[-500:] if len(result.stdout) > 500 else result.stdout,
            result.stderr[-500:] if len(result.stderr) > 500 else result.stderr,
        )
    except subprocess.TimeoutExpired:
        return (sid, -1, "", "TIMEOUT")
    except Exception as exc:
        return (sid, -1, "", str(exc))


def run_parallel(
    tasks: list[tuple[str, str, int | None]],
    *,
    label: str,
    max_workers: int,
    log_dir: Path,
    worker_fn: Callable = _run_subject_cmd,
) -> dict:
    """Run subject tasks in parallel and collect success/failed lists."""
    results = {"success": [], "failed": []}
    augmented = [(sid, cmd, timeout, str(log_dir), label.replace(" ", "_")) for sid, cmd, timeout in tasks]

    print(f"  Launching {len(tasks)} tasks with {max_workers} parallel workers...")
    print(f"  Logs -> {log_dir}/")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(worker_fn, task): task[0] for task in augmented}
        for future in as_completed(futures):
            sid = futures[future]
            try:
                sid_ret, rc, _, stderr = future.result()
                if rc == 0:
                    results["success"].append(sid_ret)
                    print(f"  ✅ {sid_ret}: done")
                else:
                    results["failed"].append(sid_ret)
                    print(f"  ❌ {sid_ret}: FAILED (rc={rc})")
                    if stderr.strip():
                        print(f"     {stderr.strip()[:200]}")
            except Exception as exc:
                results["failed"].append(sid)
                print(f"  ❌ {sid}: EXCEPTION {exc}")

    print(f"\n  {label}: ✅ {len(results['success'])} / ❌ {len(results['failed'])}")
    return results
