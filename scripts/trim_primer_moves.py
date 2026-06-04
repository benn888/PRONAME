#!/usr/bin/env python3
import argparse
import sys
from typing import Dict, Tuple
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Recalibrate the mv tag of reads after trimming with cutadapt.\n"
            "The script takes a BAM file containing mv tags (from dorado --emit-moves) "
            "and a trimmed FASTQ file, and produces a new BAM file where:\n"
            "  - the read sequence and qualities come from the trimmed FASTQ,\n"
            "  - the mv tag is truncated and reindexed to match the trimmed sequence.\n"
            "Note: the first integer of the mv tag is treated as an offset and is NOT "
            "included in the sum of moves."
        )
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM file (from dorado basecaller with --emit-moves)."
    )
    parser.add_argument(
        "--fastq",
        required=True,
        help="FASTQ file trimmed by cutadapt."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output BAM file with adjusted mv tags."
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Strict mode: if any issue is encountered (missing sequence, "
            "inconsistent mv tag, etc.), the script stops immediately."
        )
    )
    parser.add_argument(
        "--max-error-reads",
        type=int,
        default=10000,
        help=(
            "Maximum number of reads that may be skipped because move recalibration failed. "
            "If exceeded, the script exits with an error. Use -1 to disable this safeguard."
        )
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=0,
        help=(
            "Refresh the progress bar every N trimmed reads encountered. "
            "If set to 0, the interval is automatically calculated as total_trimmed_reads / 200. "
            "Use -1 to disable progress reporting. [Default: 0]"
        )
    )
    parser.add_argument(
        "--progress-steps",
        type=int,
        default=200,
        help=(
            "Number of progress bar refreshes when --progress-interval is set to 0. "
            "For example, 200 means that the progress interval will be approximately "
            "1/200 of the total number of trimmed reads, i.e. one update every 0.5%. "
            "[Default: 200]"
        )
    )
    return parser.parse_args()


def read_trimmed_fastq(path: str) -> Dict[str, Tuple[str, str]]:
    """Read a FASTQ file (one sequence line per entry) and return
    a dict {qname -> (sequence, quality_string)}.
    """
    trimmed = {}
    n = 0
    with open(path, "r") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            if not qual:
                raise RuntimeError(
                    f"Truncated or malformed FASTQ around entry {n} in {path}"
                )

            if not header.startswith("@"):
                raise RuntimeError(
                    f"Unexpected FASTQ header (does not start with '@'): {header.strip()}"
                )

            qname = header[1:].strip().split()[0]
            seq = seq.strip()
            qual = qual.strip()

            if len(seq) != len(qual):
                raise RuntimeError(
                    f"Inconsistent sequence/quality length for {qname}: "
                    f"{len(seq)} nt vs {len(qual)} scores"
                )

            trimmed[qname] = (seq, qual)
            n += 1

    print(f"[INFO] Reading trimmed FASTQ: {path}")
    print(f"[INFO] {n} trimmed reads loaded.")
    return trimmed


def phred_string_to_scores(qual: str):
    """Convert a Phred+33 quality string to a list of integers."""
    return [ord(c) - 33 for c in qual]


def print_progress_bar(
    current: int,
    total: int,
    processed: int,
    errors: int,
    width: int = 40,
    force_newline: bool = False,
):
    """Print a visual terminal progress bar on stderr."""
    if total <= 0:
        pct = 0.0
        filled = 0
    else:
        pct = current / total
        if pct > 1:
            pct = 1.0
        filled = int(width * pct)

    bar = "█" * filled + " " * (width - filled)
    percent = pct * 100

    end = "\n" if force_newline else "\r"

    print(
        f"Move-table adjustment: [{bar}] "
        f"{percent:6.2f}%  "
        f"{current}/{total} trimmed reads seen | "
        f"adjusted={processed} | errors={errors}",
        end=end,
        file=sys.stderr,
        flush=True,
    )


def adjust_mv_tag(
    original_mv,
    original_len: int,
    left_trim: int,
    trimmed_len: int,
    qname: str,
    strict: bool,
):
    """
    Adjust the mv array to match the trimmed sequence.

    IMPORTANT:
      - mv is an array of "moves" (0, 1, 2, ...) stored in a B:c tag
      - the FIRST integer is an OFFSET (e.g. 6 in mv:B:c,6,1,0,0,...)
      - the actual moves start at the 2nd integer (mv[1:])
      - sum(real_moves) == length of the original sequence
    """
    try:
        mv_list = list(original_mv)
    except TypeError:
        msg = f"[ERROR] mv tag for {qname} is not iterable: {type(original_mv)}"
        if strict:
            raise RuntimeError(msg)
        else:
            print(msg, file=sys.stderr)
            return None

    if not mv_list:
        msg = f"[ERROR] Empty mv tag for {qname}."
        if strict:
            raise RuntimeError(msg)
        else:
            print(msg, file=sys.stderr)
            return None

    # First element = offset, not counted in base sum
    offset = mv_list[0]
    moves = mv_list[1:]

    total_bases_from_moves = sum(moves)
    if total_bases_from_moves != original_len:
        msg = (
            f"[ERROR] Sum of real moves ({total_bases_from_moves}) "
            f"!= sequence length ({original_len}) for {qname} "
            f"(initial offset = {offset})."
        )
        if strict:
            raise RuntimeError(msg)
        else:
            print(msg, file=sys.stderr)
            return None

    keep_start = left_trim
    keep_end = left_trim + trimmed_len

    new_moves = []
    cum_base = 0  # number of bases consumed before the current event

    for v in moves:
        base_start = cum_base
        base_end = cum_base + v  # exclusive
        cum_base = base_end

        # Entirely before trimming window
        if base_end <= keep_start:
            continue

        # Entirely after trimming window
        if base_start >= keep_end:
            if v > 0:
                break
            else:
                continue

        # Overlaps trimming window
        if v == 0:
            # Move 0: no base advance
            if keep_start <= base_start < keep_end:
                new_moves.append(0)
            continue

        overlap_start = max(base_start, keep_start)
        overlap_end = min(base_end, keep_end)
        keep_v = overlap_end - overlap_start

        if keep_v > 0:
            new_moves.append(keep_v)

    if not new_moves:
        msg = (
            f"[ERROR] mv tag empty after trimming for {qname}. "
            f"This usually means that the retained trimmed sequence does not map to any base-emitting move interval in the original mv table."
        )
        if strict:
            raise RuntimeError(msg)
        else:
            print(msg, file=sys.stderr)
            return None

    if sum(new_moves) != trimmed_len:
        msg = (
            f"[ERROR] Sum of adjusted moves ({sum(new_moves)}) "
            f"!= trimmed length ({trimmed_len}) for {qname}."
        )
        if strict:
            raise RuntimeError(msg)
        else:
            print(msg, file=sys.stderr)
            return None

    # Rebuild final mv tag: same offset, then trimmed moves
    new_mv_full = [offset] + new_moves
    return new_mv_full


def main():
    args = parse_args()

    if args.max_error_reads < -1:
        raise RuntimeError(
            f"[ERROR] --max-error-reads must be -1 or a non-negative integer, got {args.max_error_reads}."
        )

    if args.progress_interval < -1:
        raise RuntimeError(
            f"[ERROR] --progress-interval must be -1, 0, or a positive integer, got {args.progress_interval}."
        )

    if args.progress_steps < 1:
        raise RuntimeError(
            f"[ERROR] --progress-steps must be a positive integer, got {args.progress_steps}."
        )

    trimmed_reads = read_trimmed_fastq(args.fastq)

    print(f"[INFO] Opening input BAM: {args.bam}")
    # Unaligned BAM -> no @SQ -> check_sq=False
    in_bam = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
    out_bam = pysam.AlignmentFile(args.output, "wb", header=in_bam.header)
    print(f"[INFO] Output BAM: {args.output}")

    total_reads = 0
    processed_reads = 0
    skipped_missing_in_fastq = 0
    skipped_errors = 0
    total_trimmed_reads = len(trimmed_reads)
    seen_trimmed_reads = 0
    last_progress_report = 0

    if args.progress_interval == 0:
        progress_interval = max(1, total_trimmed_reads // args.progress_steps)
    else:
        progress_interval = args.progress_interval

    for read in in_bam:
        total_reads += 1

        qname = read.query_name

        if qname not in trimmed_reads:
            # Read discarded by cutadapt -> do not keep it
            skipped_missing_in_fastq += 1
            continue

        seen_trimmed_reads += 1

        if (
            progress_interval > 0
            and (
                seen_trimmed_reads == 1
                or seen_trimmed_reads - last_progress_report >= progress_interval
                or seen_trimmed_reads == total_trimmed_reads
            )
        ):
            print_progress_bar(
                current=seen_trimmed_reads,
                total=total_trimmed_reads,
                processed=processed_reads,
                errors=skipped_errors,
            )
            last_progress_report = seen_trimmed_reads

        trimmed_seq, trimmed_qual_str = trimmed_reads[qname]

        if len(trimmed_seq) == 0:
            print(f"[WARNING] Skipping {qname}: trimmed sequence is empty.", file=sys.stderr)
            skipped_errors += 1
            continue

        if len(trimmed_seq) != len(trimmed_qual_str):
            print(
                f"[WARNING] Skipping {qname}: trimmed sequence/quality length mismatch "
                f"({len(trimmed_seq)} vs {len(trimmed_qual_str)}).",
                file=sys.stderr
            )
            skipped_errors += 1
            continue

        original_seq = read.query_sequence
        if original_seq is None:
            msg = f"[ERROR] Read {qname} has no sequence in input BAM."
            if args.strict:
                raise RuntimeError(msg)
            else:
                print(msg, file=sys.stderr)
                skipped_errors += 1
                continue

        # Locate trimmed sequence within the original sequence
        pos = original_seq.find(trimmed_seq)
        if pos == -1:
            msg = (
                f"[ERROR] Unable to locate trimmed sequence for {qname} "
                f"in the original sequence.\n"
                f"  len(original)={len(original_seq)}, len(trimmed)={len(trimmed_seq)}"
            )
            if args.strict:
                raise RuntimeError(msg)
            else:
                print(msg, file=sys.stderr)
                skipped_errors += 1
                continue

        left_trim = pos
        trimmed_len = len(trimmed_seq)
        original_len = len(original_seq)

        # Update sequence and qualities
        read.query_sequence = trimmed_seq
        read.query_qualities = phred_string_to_scores(trimmed_qual_str)

        # Retrieve original mv tag
        try:
            original_mv = read.get_tag("mv")
        except KeyError:
            msg = f"[ERROR] Read {qname} has no mv tag in input BAM."
            if args.strict:
                raise RuntimeError(msg)
            else:
                print(msg, file=sys.stderr)
                skipped_errors += 1
                continue

        # Adjust mv (preserving offset + moves logic)
        new_mv = adjust_mv_tag(
            original_mv=original_mv,
            original_len=original_len,
            left_trim=left_trim,
            trimmed_len=trimmed_len,
            qname=qname,
            strict=args.strict,
        )

        if new_mv is None:
            skipped_errors += 1
            continue

        # *** IMPORTANT ***
        # Try to preserve the original Python type of mv (often array('b')),
        # so that pysam can correctly infer the B:c type without forcing value_type.
        try:
            mv_value = type(original_mv)(new_mv)
        except Exception:
            # Fallback: simple list of ints, pysam will choose a default type
            mv_value = list(new_mv)

        # Do NOT specify value_type, let pysam infer it
        read.set_tag("mv", mv_value)

        out_bam.write(read)
        processed_reads += 1

    in_bam.close()
    out_bam.close()

    del trimmed_reads

    if progress_interval > 0:
        print_progress_bar(
            current=seen_trimmed_reads,
            total=total_trimmed_reads,
            processed=processed_reads,
            errors=skipped_errors,
            force_newline=True,
        )

    print("\n[INFO] Finished.")
    print(f"[INFO] Total reads in input BAM: {total_reads}")
    print(f"[INFO] Reads present in trimmed FASTQ and processed: {processed_reads}")
    print(f"[INFO] Reads skipped because absent from trimmed FASTQ: {skipped_missing_in_fastq}")
    print(f"[INFO] Reads skipped due to errors: {skipped_errors}")

    if processed_reads == 0:
        print(
            "[WARNING] No reads were written to the output BAM. "
            "Check that read IDs in the trimmed FASTQ and input BAM match.",
            file=sys.stderr,
        )

    if args.max_error_reads >= 0 and skipped_errors > args.max_error_reads:
        raise RuntimeError(
            f"[ERROR] Too many reads were skipped due to move-recalibration errors: "
            f"{skipped_errors} > allowed maximum ({args.max_error_reads})."
        )


if __name__ == "__main__":
    main()
