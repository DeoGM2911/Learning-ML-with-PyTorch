from celery import Celery

celery = Celery(
    'tasks',
    broker='redis://localhost:6379/0',
    result_backend='redis://localhost:6379/0'
)

@celery.task
def add(x, y):
    return x + y
