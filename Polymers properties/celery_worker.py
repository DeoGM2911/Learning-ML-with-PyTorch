# celery_worker.py
from app import create_app, celery              # `celery` is bound inside create_app

app = create_app()                              # creates Flask + Celery