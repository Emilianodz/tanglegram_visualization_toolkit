"""Logging setup: keeps the [-] / [+] / [!] / [i] message style."""
import logging
import sys


def configurar_logging(silencioso=False):
    """Send the package messages to stdout with no extra formatting."""
    logger = logging.getLogger("tanglegram_utils")
    logger.handlers.clear()
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(handler)
    logger.setLevel(logging.WARNING if silencioso else logging.INFO)
    logger.propagate = False
    return logger
