"""Run the flat-torus G0--G3 refinement audit and print JSON."""
import json
from src.flat_torus_refinement import flat_torus_refinement_audit

if __name__ == "__main__":
    report = flat_torus_refinement_audit(
        levels=(4, 8, 16),
        dimension=2,
        side_length=1.0,
        mode=(1, 0),
    )
    print(json.dumps(report, indent=2, sort_keys=True))
