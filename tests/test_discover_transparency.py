from pathlib import Path


def test_discover_has_no_silent_exception_swallowing():
    text = Path("src/phylofoundry/tasks/discover.py").read_text()
    assert "except Exception:\n            pass" not in text
