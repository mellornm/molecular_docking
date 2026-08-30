"""
Módulo de Notificações por E-mail para Etapas Demoradas do Pipeline de Docking e Dinâmica Molecular.
Utiliza smtplib (com suporte a SMTP_SSL na porta 465) e python-dotenv para envio de alertas automáticos.
"""

import os
import ssl
import smtplib
import socket
from datetime import datetime
from email.message import EmailMessage
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

try:
    from dotenv import load_dotenv

    # Carrega o .env localizado na raiz do projeto
    _PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
    _ENV_PATH = _PROJECT_ROOT / ".env"
    if _ENV_PATH.exists():
        load_dotenv(dotenv_path=_ENV_PATH)
    else:
        load_dotenv()
except ImportError:
    pass


def get_email_config() -> Dict[str, Any]:
    """
    Carrega e retorna as configurações de SMTP a partir das variáveis de ambiente (.env).
    """
    smtp_server = (
        os.getenv("SMTP_SERVER") or os.getenv("EMAIL_HOST") or "smtp.gmail.com"
    ).strip()

    port_str = os.getenv("SMTP_PORT") or os.getenv("EMAIL_PORT") or "465"
    try:
        smtp_port = int(port_str.strip())
    except ValueError:
        smtp_port = 465

    smtp_user = (
        os.getenv("SMTP_USER")
        or os.getenv("EMAIL_USER")
        or os.getenv("EMAIL_SENDER")
        or ""
    ).strip()

    smtp_password = (
        os.getenv("SMTP_PASSWORD") or os.getenv("EMAIL_PASSWORD") or ""
    ).strip()

    email_receiver = (
        os.getenv("EMAIL_RECEIVER") or os.getenv("EMAIL_TO") or smtp_user
    ).strip()

    use_ssl_env = os.getenv("SMTP_USE_SSL", "true").strip().lower()
    use_ssl = use_ssl_env in ("true", "1", "yes", "t")

    enabled_env = os.getenv("EMAIL_NOTIFICATIONS_ENABLED", "true").strip().lower()
    enabled = enabled_env in ("true", "1", "yes", "t")

    is_configured = bool(smtp_user and smtp_password)

    return {
        "smtp_server": smtp_server,
        "smtp_port": smtp_port,
        "smtp_user": smtp_user,
        "smtp_password": smtp_password,
        "email_receiver": email_receiver,
        "use_ssl": use_ssl,
        "enabled": enabled,
        "is_configured": is_configured,
    }


def is_email_configured() -> bool:
    """Verifica se as credenciais de e-mail estão configuradas no .env."""
    cfg = get_email_config()
    return cfg["is_configured"] and cfg["enabled"]


def format_duration(seconds: Optional[float]) -> str:
    """Formata a duração em segundos para um formato legível (ex: 2h 15m 30s)."""
    if seconds is None or seconds < 0:
        return "N/A"

    total_seconds = int(round(seconds))
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    secs = total_seconds % 60

    parts = []
    if hours > 0:
        parts.append(f"{hours}h")
    if minutes > 0 or hours > 0:
        parts.append(f"{minutes}m")
    parts.append(f"{secs}s")
    return " ".join(parts)


def send_email_alert(
    step_name: str,
    status: str = "success",
    duration_seconds: Optional[float] = None,
    details: Optional[Dict[str, Any]] = None,
    error_message: Optional[str] = None,
    console_logger: Optional[Any] = None,
) -> Tuple[bool, str]:
    """
    Envia um e-mail de notificação com informações detalhadas da etapa finalizada.

    :param step_name: Nome da etapa (ex: "Opção 9: Pós-processamento e MM-PBSA da DM")
    :param status: "success", "error", ou "warning"
    :param duration_seconds: Tempo de execução em segundos
    :param details: Dicionário contendo métricas e informações adicionais (ex: RMSD, Score, ΔG)
    :param error_message: Mensagem detalhada de erro caso status seja "error"
    :param console_logger: Instância de rich Console (opcional) para logs
    :return: (sucesso: bool, mensagem: str)
    """
    cfg = get_email_config()

    if not cfg["enabled"]:
        msg = "Notificações por e-mail desativadas via EMAIL_NOTIFICATIONS_ENABLED."
        if console_logger:
            console_logger.print(f"[dim]{msg}[/dim]")
        return False, msg

    if not cfg["is_configured"]:
        msg = (
            "Notificações de e-mail não configuradas no arquivo .env. "
            "Defina SMTP_USER e SMTP_PASSWORD no seu .env para receber alertas automáticos."
        )
        if console_logger:
            console_logger.print(f"[yellow]⚠ {msg}[/yellow]")
        return False, msg

    # Configuração de data e hora local
    now = datetime.now()
    now_str = now.strftime("%d/%m/%Y %H:%M:%S")
    hostname = socket.gethostname()

    # Montagem do Título / Assunto
    if status.lower() == "success":
        subject_icon = "✓"
        status_label = "Concluída com Sucesso"
        badge_color = "#28a745"
    elif status.lower() == "warning":
        subject_icon = "⚠"
        status_label = "Concluída com Alertas"
        badge_color = "#ffc107"
    else:
        subject_icon = "✗"
        status_label = "Falha na Execução"
        badge_color = "#dc3545"

    subject = (
        f"[Molecular Docking] {subject_icon} Etapa Pronta: {step_name} [{status_label}]"
    )

    duration_str = format_duration(duration_seconds)

    # Monta texto simples (Plain Text Fallback)
    lines = [
        f"=== ALERTA DE EXECUÇÃO: {step_name.upper()} ===",
        f"Status: {status_label}",
        f"Data/Hora: {now_str}",
        f"Duração: {duration_str}",
        f"Servidor/Máquina: {hostname}",
        "",
    ]

    if details:
        lines.append("--- DETALHES & RESULTADOS ---")
        for k, v in details.items():
            lines.append(f"• {k}: {v}")
        lines.append("")

    if error_message:
        lines.append("--- DETALHES DO ERRO ---")
        lines.append(str(error_message))
        lines.append("")

    lines.append("Mensagem gerada automaticamente pelo Pipeline de Docking Molecular.")
    plain_content = "\n".join(lines)

    # Monta HTML estilizado
    details_html = ""
    if details:
        details_rows = "".join(
            f"<tr><td style='padding: 8px 12px; border-bottom: 1px solid #e9ecef; font-weight: 600; color: #495057;'>{k}</td>"
            f"<td style='padding: 8px 12px; border-bottom: 1px solid #e9ecef; color: #212529; font-family: monospace;'>{v}</td></tr>"
            for k, v in details.items()
        )
        details_html = f"""
        <div style="margin-top: 20px;">
            <h3 style="color: #343a40; font-size: 16px; margin-bottom: 10px; border-bottom: 2px solid #dee2e6; padding-bottom: 5px;">
                📊 Resultados e Métricas Principais
            </h3>
            <table style="width: 100%; border-collapse: collapse; background-color: #f8f9fa; border-radius: 6px; overflow: hidden;">
                <tbody>
                    {details_rows}
                </tbody>
            </table>
        </div>
        """

    error_html = ""
    if error_message:
        error_html = f"""
        <div style="margin-top: 20px; background-color: #fff5f5; border: 1px solid #f5c6cb; border-radius: 6px; padding: 15px;">
            <h3 style="color: #721c24; font-size: 15px; margin-top: 0; margin-bottom: 8px;">
                ⚠️ Detalhes do Erro
            </h3>
            <pre style="margin: 0; font-family: monospace; font-size: 13px; color: #721c24; white-space: pre-wrap; word-break: break-word;">{error_message}</pre>
        </div>
        """

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>{subject}</title>
    </head>
    <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; background-color: #f4f6f9; margin: 0; padding: 20px; color: #333;">
        <table align="center" border="0" cellpadding="0" cellspacing="0" width="600" style="background-color: #ffffff; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); overflow: hidden; border: 1px solid #e1e4e8;">
            <!-- Header -->
            <tr>
                <td style="background-color: #1a365d; padding: 24px; text-align: center; color: #ffffff;">
                    <h1 style="margin: 0; font-size: 20px; letter-spacing: 0.5px;">🧬 Pipeline de Docking Molecular</h1>
                    <p style="margin: 6px 0 0 0; font-size: 13px; color: #cbd5e0;">Sistema Automatizado de Notificações</p>
                </td>
            </tr>
            <!-- Main Content -->
            <tr>
                <td style="padding: 24px;">
                    <div style="display: flex; align-items: center; justify-content: space-between; margin-bottom: 18px; border-bottom: 1px solid #edf2f7; padding-bottom: 14px;">
                        <span style="display: inline-block; background-color: {badge_color}; color: #ffffff; padding: 6px 14px; font-size: 12px; font-weight: bold; border-radius: 20px; text-transform: uppercase;">
                            {status_label}
                        </span>
                        <span style="font-size: 13px; color: #718096; float: right;">
                            ⏱ <strong>{duration_str}</strong> de execução
                        </span>
                    </div>

                    <h2 style="font-size: 18px; color: #2d3748; margin-top: 0; margin-bottom: 12px;">
                        Etapa: <span style="color: #2b6cb0;">{step_name}</span>
                    </h2>
                    
                    <p style="font-size: 14px; color: #4a5568; line-height: 1.5; margin: 0 0 16px 0;">
                        A execução desta etapa foi finalizada em <strong>{now_str}</strong> na estação <code>{hostname}</code>.
                    </p>

                    {details_html}
                    {error_html}

                    <div style="margin-top: 25px; padding-top: 15px; border-top: 1px solid #edf2f7; font-size: 12px; color: #a0aec0; text-align: center;">
                        Notificação enviada para <code>{cfg["email_receiver"]}</code> via SMTP SSL ({cfg["smtp_server"]}:{cfg["smtp_port"]}).
                    </div>
                </td>
            </tr>
        </table>
    </body>
    </html>
    """

    # Construção do objeto EmailMessage
    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"] = cfg["smtp_user"]
    msg["To"] = cfg["email_receiver"]
    msg.set_content(plain_content)
    msg.add_alternative(html_content, subtype="html")

    try:
        if cfg["smtp_port"] == 465 or cfg["use_ssl"]:
            context = ssl.create_default_context()
            with smtplib.SMTP_SSL(
                host=cfg["smtp_server"],
                port=cfg["smtp_port"],
                context=context,
                timeout=30,
            ) as server:
                server.login(cfg["smtp_user"], cfg["smtp_password"])
                server.send_message(msg)
        else:
            # Fallback para portas como 587 com STARTTLS
            with smtplib.SMTP(
                host=cfg["smtp_server"],
                port=cfg["smtp_port"],
                timeout=30,
            ) as server:
                context = ssl.create_default_context()
                server.starttls(context=context)
                server.login(cfg["smtp_user"], cfg["smtp_password"])
                server.send_message(msg)

        success_msg = (
            f"E-mail de alerta enviado com sucesso para {cfg['email_receiver']}!"
        )
        if console_logger:
            console_logger.print(f"[bold green]📧 {success_msg}[/bold green]")
        return True, success_msg

    except Exception as e:
        err_msg = f"Falha ao enviar e-mail de alerta ({type(e).__name__}): {e}"
        if console_logger:
            console_logger.print(
                f"[bold yellow]⚠ Aviso de Notificação:[/bold yellow] {err_msg}"
            )
        return False, err_msg


def test_email_connection(console_logger: Optional[Any] = None) -> Tuple[bool, str]:
    """
    Envia um e-mail de teste para verificar se as credenciais e o servidor SMTP estão corretos.
    """
    cfg = get_email_config()
    if not cfg["is_configured"]:
        msg = "Credenciais ausentes no .env. Configure SMTP_USER e SMTP_PASSWORD."
        if console_logger:
            console_logger.print(f"[bold red]✗ {msg}[/bold red]")
        return False, msg

    test_details = {
        "Servidor SMTP": f"{cfg['smtp_server']}:{cfg['smtp_port']}",
        "Remetente": cfg["smtp_user"],
        "Destinatário": cfg["email_receiver"],
        "Modo de Segurança": "SMTP_SSL (Porta 465)"
        if (cfg["use_ssl"] or cfg["smtp_port"] == 465)
        else f"STARTTLS (Porta {cfg['smtp_port']})",
        "Status do Teste": "Configuração validada com sucesso!",
    }

    return send_email_alert(
        step_name="Teste de Conectividade de E-mail",
        status="success",
        duration_seconds=1.2,
        details=test_details,
        console_logger=console_logger,
    )
