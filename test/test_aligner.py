import subprocess
import sys
from pathlib import Path


def run_cli(fasta_name: str, algorithm: str | None = None):
    """
    Runs: python3 src/PairwiseAlignment/align.py <fasta> [algorithm]
    Returns CompletedProcess
    """
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "src" / "PairwiseAlignment" / "align.py"
    fasta = Path(__file__).parent / "testdata" / fasta_name

    cmd = [sys.executable, str(script), str(fasta)]
    if algorithm:
        cmd.append(algorithm)

    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=repo_root,
    )


def test_identical_default_nw_succeeds():
    res = run_cli("identical.fasta")

    assert res.returncode == 0, f"STDERR:\n{res.stderr}\nSTDOUT:\n{res.stdout}"
    assert "Running Needleman-Wunsch (Global Alignment)" in res.stdout
    assert "Final Alignment Score:" in res.stdout
    assert "S1: ACGTACGT" in res.stdout
    assert "S2: ACGTACGT" in res.stdout


def test_identical_sw_succeeds():
    res = run_cli("identical.fasta", "sw")

    assert res.returncode == 0, f"STDERR:\n{res.stderr}\nSTDOUT:\n{res.stdout}"
    assert "Running Smith-Waterman (Local Alignment)" in res.stdout
    assert "Final Alignment Score: 16" in res.stdout


def test_mismatch_default_nw_score():
    res = run_cli("mismatch.fasta")

    assert res.returncode == 0, f"STDERR:\n{res.stderr}\nSTDOUT:\n{res.stdout}"
    assert "Running Needleman-Wunsch (Global Alignment)" in res.stdout
    assert "Final Alignment Score: 10" in res.stdout
    assert "S1: ACGTACGT" in res.stdout
    assert "S2: ACGTTCGA" in res.stdout


def test_edge_empty_like_reports_error_message():
    res = run_cli("edge_empty_like.fasta")

    combined = f"{res.stdout}\n{res.stderr}"
    assert "Error: Cannot align. One or both sequences are empty or contain invalid characters." in combined
