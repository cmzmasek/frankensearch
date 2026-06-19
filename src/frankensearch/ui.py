"""Console output helpers (Rich) for friendly, consistent messaging."""

from __future__ import annotations

from rich.console import Console
from rich.panel import Panel

console = Console()
err_console = Console(stderr=True)


def info(message: str) -> None:
    console.print(message)


def warn(message: str) -> None:
    err_console.print(f"[bold yellow]Warning:[/] {message}")


def error_panel(title: str, message: str, hint: str | None = None) -> None:
    body = message
    if hint:
        body += f"\n\n[bold]Try:[/] {hint}"
    err_console.print(
        Panel(body, title=f"[bold red]{title}[/]", border_style="red", expand=False)
    )
