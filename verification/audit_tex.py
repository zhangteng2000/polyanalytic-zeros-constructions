"""Check labels, equation references, citations, and the agreed TeX conventions."""
from pathlib import Path
import collections
import json
import re
import sys

def audit(path):
    source=path.read_text(encoding="utf-8")
    # Strip unescaped comments.
    source=re.sub(r"(?<!\\)%[^\n]*","",source)
    labels=re.findall(r"\\label\{([^}]+)\}",source)
    refs=re.findall(r"\\(?:eqref|ref|pageref|autoref)\{([^}]+)\}",source)
    citations=[k.strip() for group in re.findall(r"\\cite(?:\[[^\]]*\]){0,2}\{([^}]+)\}",source)
               for k in group.split(",")]
    bib=re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}",source)
    equation_labels=[]
    for block in re.finditer(r"\\begin\{(equation|align|gather|multline)\}([\s\S]*?)\\end\{\1\}",source):
        equation_labels.extend(re.findall(r"\\label\{([^}]+)\}",block.group(2)))
    optional=[]
    for item in re.finditer(r"\\begin\{([^}]+)\}\[([^\n]*)",source):
        if item.group(1)=="lem" and item.group(2).startswith(r"{\cite["):continue
        if item.group(1) in ("figure","tikzpicture","scope"):continue
        optional.append(item.group(0))
    report={
        "file":path.name,"label_count":len(labels),"numbered_equations":len(equation_labels),
        "bibliography_count":len(bib),
        "duplicate_labels":[k for k,v in collections.Counter(labels).items() if v>1],
        "undefined_references":sorted(set(refs)-set(labels)),
        "uncited_numbered_equations":sorted(set(equation_labels)-set(refs)),
        "undefined_citations":sorted(set(citations)-set(bib)),
        "uncited_bibliography_items":sorted(set(bib)-set(citations)),
        "disallowed_optional_environment_titles":optional,
        "proofs_containing_external_citations":len([
            b for b in re.findall(r"\\begin\{proof\}([\s\S]*?)\\end\{proof\}",source)
            if r"\cite" in b]),
        "non_ascii_characters":sorted(set(c for c in source if ord(c)>127)),
    }
    keys=("duplicate_labels","undefined_references","uncited_numbered_equations",
          "undefined_citations","uncited_bibliography_items","disallowed_optional_environment_titles",
          "proofs_containing_external_citations")
    report["status"]="PASS" if not any(report[k] for k in keys) else "NEEDS_REVIEW"
    return report

if __name__=="__main__":
    paths=[Path(x) for x in sys.argv[1:]]
    if not paths:paths=sorted((Path(__file__).resolve().parents[1]/"paper").glob("*.tex"))
    reports=[audit(path) for path in paths]
    print(json.dumps(reports,indent=2,ensure_ascii=True))
    (Path(__file__).resolve().parent/"tex_audit.json").write_text(
        json.dumps(reports,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    sys.exit(0 if all(r["status"]=="PASS" for r in reports) else 1)
