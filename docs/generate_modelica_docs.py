#!/usr/bin/env python3

import re
from pathlib import Path
import shutil


PACKAGE_ROOT = Path("../src/rtctools_channel_flow/modelica/Deltares")
OUTPUT_DIR = Path("source/generated")




if OUTPUT_DIR.exists():
    shutil.rmtree(OUTPUT_DIR)

OUTPUT_DIR.mkdir(parents=True)


OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


MATH_RE = re.compile(r"<math>(.*?)</math>", re.DOTALL | re.IGNORECASE)
INLINE_CODE_RE = re.compile(
    r"<code>([^<>\n]+)</code>",
    re.IGNORECASE,
)

BLOCK_CODE_RE = re.compile(
    r"<code>(.*?)</code>",
    re.DOTALL | re.IGNORECASE,
)
MODEL_RE = re.compile(
    r"\b(model|package|block|record)\s+([A-Za-z0-9_]+)"
)

DOC_RE = re.compile(
    r"Documentation\s*\(\s*info\s*=\s*\"(.*?)\"\s*\)",
    re.DOTALL,
)


def clean_html_fragment(text):
    """Convert a small subset of Modelica HTML documentation to RST."""

    # Remove wrapper tags
    text = re.sub(r"</?html>", "", text, flags=re.IGNORECASE)

    # Paragraphs and line breaks
    text = re.sub(r"<p>", "", text, flags=re.IGNORECASE)
    text = re.sub(r"</p>", "\n\n", text, flags=re.IGNORECASE)
    text = re.sub(r"<br\s*/?>", "\n", text, flags=re.IGNORECASE)

    # Lists
    text = re.sub(r"<ul>", "\n", text, flags=re.IGNORECASE)
    text = re.sub(r"</ul>", "\n", text, flags=re.IGNORECASE)
    text = re.sub(r"<li>", "- ", text, flags=re.IGNORECASE)
    text = re.sub(r"</li>", "\n", text, flags=re.IGNORECASE)
    text = INLINE_CODE_RE.sub(r"``\1``", text)

    # Inline formatting
    text = re.sub(
        r"<em>(.*?)</em>",
        r"*\1*",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    text = re.sub(
        r"<i>(.*?)</i>",
        r"*\1*",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    text = re.sub(
        r"<b>(.*?)</b>",
        r"**\1**",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    text = re.sub(
        r"<strong>(.*?)</strong>",
        r"**\1**",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )

    text = re.sub(
        r"<code>([^<>\n]+)</code>",
        r"``\1``",
        text,
        flags=re.IGNORECASE,
    )
    '''
    # Inline code
    text = re.sub(
        r"<code>(.*?)</code>",
        r"``\1``",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    '''




    return text.strip()


def html_to_rst(doc):
    """Convert Modelica Documentation(info=...) HTML to RST.

    Supports:
    - <p>...</p>
    - <em>...</em>
    - <ul><li>...</li></ul>
    - <math>...</math>
    - <code>...</code>
    """

    rst = []
    pos = 0

    # Process both math and code blocks in document order
    '''
    block_re = re.compile(
        r"(<math>.*?</math>|<code>.*?</code>)",
        re.DOTALL | re.IGNORECASE,
    )
    '''

    block_re = re.compile(
        r"(<math>.*?</math>|<code>.*?\n.*?</code>)",
        re.DOTALL | re.IGNORECASE,
        )


    for match in block_re.finditer(doc):
        before = doc[pos:match.start()]
        before = clean_html_fragment(before)

        if before:
            rst.append(before)
            rst.append("")

        block = match.group(1)

        math_match = MATH_RE.fullmatch(block.strip())
        code_match = BLOCK_CODE_RE.fullmatch(block.strip())

        if math_match:
            equation = math_match.group(1).strip()

            # Convert escaped Modelica backslashes to LaTeX backslashes
            equation = equation.replace("\\\\", "\\")

            rst.append(".. math::")
            rst.append("")

            for line in equation.splitlines():
                rst.append(f"   {line.strip()}")

            rst.append("")

        elif code_match:
            code = code_match.group(1).strip()

            rst.append("Example:")
            rst.append("")
            rst.append("::")
            rst.append("")

            for line in code.splitlines():
                rst.append(f"   {line.rstrip()}")

            rst.append("")

        pos = match.end()

    remaining = doc[pos:]
    remaining = clean_html_fragment(remaining)

    if remaining:
        rst.append(remaining)

    result = "\n".join(rst)

    # Remove excessive blank lines
    result = re.sub(r"\n{3,}", "\n\n", result)

    return result.strip()


def extract_documentation(text):
    match = DOC_RE.search(text)
    if not match:
        return None

    doc = match.group(1)

    # Convert Modelica string escapes
    doc = doc.replace(r"\"", '"')
    doc = doc.replace(r"\n", "\n")

    return doc


def extract_name(text, fallback):
    match = MODEL_RE.search(text)
    return match.group(2) if match else fallback



def generate_root_index(package_contents, output_dir):
    """
    Create generated/index.rst listing the top-level packages.
    """

    top_packages = set()

    for package_path in package_contents:
        if package_path:
            top_packages.add(package_path[0])

    lines = [
        "Model Documentation",
        "===================",
        "",
        ".. toctree::",
        "   :maxdepth: 2",
        "",
    ]

    for package in sorted(top_packages):
        lines.append(f"   {package}/index")

    (output_dir / "index.rst").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )



def generate_package_indexes(package_contents, output_dir):
    """
    Generate index.rst files for all package levels.

    Example generated structure:

    generated/
    ├── index.rst
    ├── Hydraulic/
    │   ├── index.rst
    │   └── Reservoir/
    │       ├── index.rst
    │       └── HomotopicPower.rst
    └── SimpleRouting/
        ├── index.rst
        └── Nodes/
            ├── index.rst
            └── Node.rst
    """

    child_packages = defaultdict(set)

    # Build parent -> child package mapping
    for package_path in package_contents:
        for i in range(len(package_path)):
            parent = package_path[:i]
            child = package_path[i]
            child_packages[parent].add(child)

    all_package_paths = set(package_contents.keys())

    for package_path in package_contents:
        for i in range(1, len(package_path)):
            all_package_paths.add(package_path[:i])

    for package_path in sorted(all_package_paths):
        package_dir = output_dir.joinpath(*package_path)
        package_dir.mkdir(parents=True, exist_ok=True)

        title = package_path[-1] if package_path else "Model Documentation"

        lines = [
            title,
            "=" * len(title),
            "",
            ".. toctree::",
            "   :maxdepth: 1",
            "",
        ]

        # Add child packages first
        for child in sorted(child_packages.get(package_path, [])):
            lines.append(f"   {child}/index")

        # Add model pages in this package
        for page in sorted(package_contents.get(package_path, [])):
            lines.append(f"   {page}")

        (package_dir / "index.rst").write_text(
            "\n".join(lines) + "\n",
            encoding="utf-8",
        )


'''
def generate_package_indexes(package_contents, output_dir):

    for package_path, pages in package_contents.items():

        package_dir = output_dir.joinpath(*package_path)
        package_dir.mkdir(parents=True, exist_ok=True)

        title = package_path[-1]

        lines = [
            title,
            "=" * len(title),
            "",
            ".. toctree::",
            "   :maxdepth: 1",
            "",
        ]

        for page in sorted(pages):
            lines.append(f"   {page}")

        (package_dir / "index.rst").write_text(
            "\n".join(lines),
            encoding="utf-8",
        )
'''

print(PACKAGE_ROOT.resolve())
print(PACKAGE_ROOT.exists())
from collections import defaultdict

package_contents = defaultdict(list)

for mo_file in PACKAGE_ROOT.rglob("*.mo"):
    text = mo_file.read_text(encoding="utf-8")

    doc = extract_documentation(text)

    if not doc:
        continue

    print(f"Found documentation in {mo_file}")
    name = extract_name(text, mo_file.stem)
    body = html_to_rst(doc).strip()

    relative_path = mo_file.relative_to(PACKAGE_ROOT)

    # Example:
    # Hydraulic/Reservoir/HomotopicPower.mo
    # becomes ("Hydraulic", "Reservoir")
    package_parts = relative_path.parts[:-1]

    package_contents[package_parts].append(name)

    rst = (
        f"{name}\n"
        f"{'=' * len(name)}\n\n"
        f"{body}\n"
    )

    target_dir = OUTPUT_DIR.joinpath(*package_parts)
    target_dir.mkdir(parents=True, exist_ok=True)

    print(f"Writing {target_dir / f'{name}.rst'}")

    (target_dir / f"{name}.rst").write_text(
        rst,
        encoding="utf-8",
    )


generate_package_indexes(package_contents, OUTPUT_DIR)
generate_root_index(package_contents, OUTPUT_DIR)

print("Generated documentation pages.")

