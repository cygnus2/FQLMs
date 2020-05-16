""" ----------------------------------------------------------------------------

    lr_notify.py - LR, May 2020

    Minimal wrapper for the pushover package to send messages to Lukas when
    runs are done.

---------------------------------------------------------------------------- """
from pushover import Client


def _send(msg, device):
    client = Client("uwxqcfqmy67k1iyfvo3t74ej9no47r", api_token="ahyydd79r5ymfc2nj91a4bcnyg9qmg")
    client.send_message(
        msg['text'],
        title=msg['title'],
        device=device
    )

def _timestamp():
    return datetime.now().strftime("%b %d %Y %H:%M:%S")

def push_finished_run(run, device='oneplus_t3_luki'):
    msg = {
        'title' : 'Run completed!',
        'text' : f"The script with {run['n_tasks']} tasks @{run['host']} finished at {_timestamp()}"
    }
    _send(msg, device=device)


def push_message(text, device='oneplus_t3_luki'):
    msg = {
        'title' : 'Run completed!',
        'text' : text
    }
    _send(msg, device=device)
