from flask import Flask
from celery import Celery

celery = Celery(__name__)  # global instance (can be imported anywhere)

def create_app():
    app = Flask(__name__)

    # --- App Configuration
    app.config.from_mapping(
        SECRET_KEY='your-secret-key',
        CELERY_BROKER_URL='redis://localhost:6379/0',
        result_backend='redis://localhost:6379/0'
    )

    # --- Init Celery with app context
    init_celery(app)

    return app

def init_celery(app=None):
    app = app or create_app()
    celery.conf.update(app.config)

    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
