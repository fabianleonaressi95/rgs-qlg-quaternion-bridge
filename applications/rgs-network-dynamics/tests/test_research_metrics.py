import numpy as np
import pytest

from competence_model import concentration_by_skill, effective_agent_count, qualified_counts
from knowledge_dynamics import update_knowledge
from resilience import retention_ratio, task_coverage, targeted_agent_order


def test_competence_concentration_and_effective_count():
    c = np.array([[1.0, 0.0], [1.0, 0.5], [0.0, 0.5]])
    np.testing.assert_allclose(qualified_counts(c, 0.5), [2, 2])
    np.testing.assert_allclose(concentration_by_skill(c), [0.5, 0.5])
    np.testing.assert_allclose(effective_agent_count(c), [2.0, 2.0])


def test_uncovered_skill_is_nan_not_zero():
    c = np.zeros((3, 2))
    assert np.isnan(concentration_by_skill(c)).all()
    assert np.isnan(effective_agent_count(c)).all()


def test_knowledge_update_preserves_bounds():
    x = np.array([[0.2], [0.8]])
    p = np.array([[0.0, 1.0], [1.0, 0.0]])
    y = update_knowledge(x, p, learning_input=0.05, eta=0.2, decay=0.01)
    assert np.all((y >= 0.0) & (y <= 1.0))


def test_knowledge_update_rejects_non_stochastic_matrix():
    with pytest.raises(ValueError):
        update_knowledge(np.ones((2, 1)) * 0.5, np.ones((2, 2)))


def test_task_coverage_and_agent_removal():
    c = np.array([[1.0, 0.0], [0.0, 1.0], [0.8, 0.8]])
    tasks = np.array([[1, 0], [0, 1], [1, 1]], dtype=bool)
    assert task_coverage(c, tasks, threshold=0.75) == 1.0
    remaining = np.array([True, False, False])
    assert task_coverage(c, tasks, threshold=0.75, available_agents=remaining) == pytest.approx(1 / 3)


def test_retention_ratio_and_unique_skill_order():
    assert retention_ratio(1.0, 0.5) == 0.5
    assert np.isnan(retention_ratio(0.0, 0.0))
    c = np.array([[1.0, 0.0], [0.0, 1.0], [0.8, 0.8]])
    assert targeted_agent_order(c, 0.75)[0] in (0, 1)
