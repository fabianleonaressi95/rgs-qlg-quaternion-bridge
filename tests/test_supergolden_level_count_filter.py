from decimal import Decimal
from scripts.filter_supergolden_level_counts import analyse, supergolden

def test_supergolden_root():
    s = supergolden()
    assert abs(s**3-s**2-1) < Decimal("1e-50")

def test_exact_vector_recurrence_passes():
    rows = [
        {"n":0,"A":1,"B":0,"C":0,"total":1},
        {"n":1,"A":1,"B":1,"C":0,"total":2},
        {"n":2,"A":2,"B":1,"C":1,"total":4},
        {"n":3,"A":3,"B":1,"C":1,"total":5},
        {"n":4,"A":4,"B":2,"C":1,"total":7},
    ]
    result = analyse(rows)
    assert result["stage1_exact_recurrence_gate_passed"]
    assert result["vector_recurrence"]["passed"]

def test_bad_vector_recurrence_fails():
    rows = [
        {"n":0,"A":1,"B":0,"C":0,"total":1},
        {"n":1,"A":1,"B":1,"C":0,"total":2},
        {"n":2,"A":2,"B":1,"C":1,"total":4},
        {"n":3,"A":3,"B":2,"C":2,"total":7},
    ]
    assert not analyse(rows)["stage1_exact_recurrence_gate_passed"]

def test_perron_frequencies_sum_to_one():
    s = supergolden()
    assert abs(s**-2+s**-3+s**-4-1) < Decimal("1e-50")
