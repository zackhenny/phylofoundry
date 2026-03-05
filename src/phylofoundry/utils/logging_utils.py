import logging
import os
import time
from contextlib import contextmanager


def setup_logging(outdir: str, level: str = "INFO") -> logging.Logger:
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    logfile = os.path.join(outdir, "logs", "phylofoundry.log")
    logger = logging.getLogger("phylofoundry")
    logger.handlers.clear()
    logger.setLevel(getattr(logging, str(level).upper(), logging.INFO))

    fmt = logging.Formatter("%(asctime)s %(levelname)s %(name)s - %(message)s")
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    fh = logging.FileHandler(logfile)
    fh.setFormatter(fmt)
    logger.addHandler(sh)
    logger.addHandler(fh)
    return logger


@contextmanager
def step_timer(logger: logging.Logger, step: str, extra: str = ""):
    t0 = time.time()
    logger.info("[step:start] %s %s", step, extra)
    try:
        yield
    finally:
        dt = time.time() - t0
        logger.info("[step:end] %s duration_sec=%.2f", step, dt)
