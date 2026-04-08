#!/usr/bin/env python3
"""Convert R or Python scripts to MyST Markdown notebooks.

Rules:
  - Lines of 6+ consecutive '#' become a cell separator: +++ (first one is dropped).
  - Consecutive lines starting with '### ' (at column 0) become markdown cells.
    If the text matches a heading pattern (Fig. 1d:, Supplementary Figure 7.,
    Supplementary Fig. 2a:, Step-3:, Figure 1., etc.) the label becomes an
    ## heading, followed by the rest as prose.
  - Everything else is gathered into ```{code-cell} <lang>``` blocks.
  - Trailing whitespace is stripped from all output lines.

Usage:
    python bin/r2myst.py path/to/script.R
    python bin/r2myst.py path/to/script.py
    python bin/r2myst.py path/to/*.R          # multiple files

Output is written next to the original file with a .md extension.
"""

import re
import sys
from pathlib import Path

SEPARATOR_RE = re.compile(r"^#{6,}\s*$")
COMMENT_RE = re.compile(r"^### (.*)$")
# Patterns promoted to ## headings:
#   "Fig. 1b: ...", "Supplementary Fig. 2a: ...",
#   "Supplementary Figure 7. ...", "Step-3: ...", "Figure 1. ..."
HEADING_RE = re.compile(
    r"^((?:Supplementary\s+)?Fig(?:ure)?\s*\.?\s*\d+[a-zA-Z]?|Step-?\d+)"
    r"[.:]\s*(.*)$"
)

# Match ggsave("filename.ext", plot_expr, ...) with a string-literal filename.
# Captures the full plot expression (e.g. "p", "p1 + p2", "(p1 | p2) / (p3 | p4)")
# by grabbing everything up to a known named argument or closing paren.
GGSAVE_RE = re.compile(
    r'^ggsave\("([^"]+)\.\w+",\s*(.+?)(?:,\s*(?:width|height|dpi|units|limitsize|scale|device|path)\s*=.*)?[,)]*$'
)
# Match saveWidget(varname, "filename.ext")
SAVEWIDGET_RE = re.compile(
    r'^saveWidget\((\w+),\s*"([^"]+)\.\w+"'
)
# Match pdf("filename.ext" ...) — string-literal form
PDF_OPEN_RE = re.compile(r'^pdf\("([^"]+)\.\w+"')
# Match any pdf(...) call (including paste0 form)
PDF_OPEN_ANY_RE = re.compile(r'^pdf\(')
# Match dev.off()
DEV_OFF_RE = re.compile(r'^dev\.off\(\)')

LANG_MAP = {
    ".R": "r",
    ".r": "r",
    ".py": "python",
}

KERNEL_MAP = {
    "r": {"name": "ir", "display_name": "R"},
    "python": {"name": "python3", "display_name": "Python 3"},
}


def detect_lang(path: Path) -> str:
    return LANG_MAP.get(path.suffix, "r")


def classify_line(line: str):
    """Return (kind, payload) for a source line.

    kind is one of: 'separator', 'comment', 'code', 'blank'
    """
    stripped = line.rstrip("\n\r")
    if SEPARATOR_RE.match(stripped):
        return ("separator", "")
    m = COMMENT_RE.match(stripped)
    if m:
        return ("comment", m.group(1))
    if stripped == "":
        return ("blank", "")
    return ("code", stripped)


def rstrip_lines(lines: list[str]) -> list[str]:
    """Strip trailing whitespace from every line."""
    return [l.rstrip() for l in lines]


def flush_code(lines: list[str], lang: str) -> list[str]:
    """Wrap accumulated code lines in a code-cell directive."""
    lines = rstrip_lines(lines)
    # Trim leading/trailing blank lines
    while lines and lines[0] == "":
        lines.pop(0)
    while lines and lines[-1] == "":
        lines.pop()
    if not lines:
        return []

    # Comment out pdf()/dev.off() and extract label from pdf("literal.ext")
    label = ""
    for i, line in enumerate(lines):
        if PDF_OPEN_ANY_RE.match(line):
            m = PDF_OPEN_RE.match(line)
            if m and not label:
                label = m.group(1)
            lines[i] = "# " + line
        elif DEV_OFF_RE.match(line):
            lines[i] = "# " + line

    # Detect trailing ggsave/saveWidget and append display variable
    # (check last non-blank, non-comment line)
    last = lines[-1]
    m = GGSAVE_RE.match(last)
    if m and not re.match(r'\w+\s*=', m.group(2)):
        if not label:
            label = m.group(1)
        lines.append(m.group(2).strip())
    else:
        m = SAVEWIDGET_RE.match(last)
        if m:
            lines.append(m.group(1))
            if not label:
                label = m.group(2)

    if label:
        opener = f"```{{code-cell #{label}}} {lang}"
    else:
        opener = f"```{{code-cell}} {lang}"
    out = [opener, *lines, "```", ""]
    return out


def format_comment(text: str) -> str:
    """Turn a comment payload into markdown, promoting figure/step labels to headings."""
    text = text.rstrip()
    m = HEADING_RE.match(text)
    if m:
        label = m.group(1)
        rest = m.group(2).strip()
        parts = [f"## {label}", ""]
        if rest:
            parts.append(rest)
            parts.append("")
        return "\n".join(parts)
    return text


def make_frontmatter(title: str, lang: str) -> list[str]:
    """Build a YAML frontmatter block with title and kernelspec."""
    kern = KERNEL_MAP[lang]
    return [
        "---",
        f"title: {title}",
        "kernelspec:",
        f"  name: {kern['name']}",
        f"  display_name: {kern['display_name']}",
        "---",
        "",
    ]


def extract_title(comment_buf: list[str]) -> tuple[str, list[str]]:
    """Pull the title from the first comment block.

    The first heading (matched by HEADING_RE) becomes the title string.
    If no heading is found, the first non-empty line is used instead.
    Returns (title, remaining_lines) where remaining_lines are the
    leftover comments that should still be emitted as markdown.
    """
    # format_comment may have produced multi-line strings; flatten
    flat: list[str] = []
    for entry in comment_buf:
        flat.extend(entry.split("\n"))

    # Look for first ## heading produced by format_comment
    title = ""
    rest_start = 0
    for i, line in enumerate(flat):
        if line.startswith("## "):
            title = line[3:].strip()
            rest_start = i + 1
            # skip blank line after heading
            if rest_start < len(flat) and flat[rest_start].strip() == "":
                rest_start += 1
            break
        if not title and line.strip():
            title = line.strip()
            rest_start = i + 1
            break

    remaining = flat[rest_start:]
    # trim leading blank lines from remaining
    while remaining and remaining[0].strip() == "":
        remaining.pop(0)
    return title, remaining


def convert(src: Path) -> Path:
    lang = detect_lang(src)
    lines = src.read_text().splitlines()

    output: list[str] = []
    code_buf: list[str] = []
    comment_buf: list[str] = []
    separator_count = 0  # drop the first +++
    frontmatter_done = False

    def drain_code():
        nonlocal code_buf
        flushed = flush_code(code_buf, lang)
        if flushed:
            output.extend(flushed)
        code_buf = []

    def drain_comments():
        nonlocal comment_buf, frontmatter_done
        if not comment_buf:
            return
        # First comment block after first separator becomes frontmatter
        if not frontmatter_done and separator_count >= 1:
            title, remaining = extract_title(comment_buf)
            if title:
                output.extend(make_frontmatter(title, lang))
                if remaining:
                    output.extend(remaining)
                    output.append("")
                frontmatter_done = True
                comment_buf = []
                return
        output.extend(comment_buf)
        output.append("")
        comment_buf = []

    for raw_line in lines:
        kind, payload = classify_line(raw_line)

        if kind == "separator":
            drain_code()
            drain_comments()
            separator_count += 1
            if separator_count > 1:
                output.append("+++")
                output.append("")
            continue

        if kind == "comment":
            # If we had code building up, flush it first
            drain_code()
            comment_buf.append(format_comment(payload))
            continue

        # 'code' or 'blank' -- belongs to a code cell
        # flush any pending comments first
        drain_comments()
        code_buf.append(raw_line.rstrip("\n\r"))

    # Flush remaining buffers
    drain_code()
    drain_comments()

    # Remove trailing blank lines
    while output and output[-1].strip() == "":
        output.pop()

    dest = src.with_suffix(".md")
    dest.write_text("\n".join(output) + "\n")
    return dest


def main():
    if len(sys.argv) < 2:
        print(__doc__.strip())
        sys.exit(1)

    paths = [Path(a) for a in sys.argv[1:]]
    for p in paths:
        if not p.exists():
            print(f"ERROR: {p} does not exist", file=sys.stderr)
            continue
        if p.suffix not in LANG_MAP:
            print(f"SKIP:  {p} (unsupported extension)", file=sys.stderr)
            continue
        dest = convert(p)
        print(f"  {p}  ->  {dest}")


if __name__ == "__main__":
    main()
