from app.tasks import celery

# No need to define anything else — just point Celery to this module
# To run:
# celery -A celery_worker.celery worker --loglevel=info