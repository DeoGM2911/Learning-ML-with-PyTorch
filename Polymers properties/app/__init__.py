# app/__init__.py
from flask import Flask
from celery import Celery

# We create an *empty* placeholder so that other modules can import it
celery = Celery(__name__, broker="redis://localhost:6379/0", backend="redis://localhost:6379/0")

def make_celery(app: Flask) -> Celery:
    """Factory that builds a Celery instance bound to `app`."""
    celery_app = Celery(
        app.import_name,
        broker=app.config["broker_url"],
        backend=app.config["result_backend"],
        include=["app.tasks"],          # tasks are auto-discovered
    )

    celery_app.conf.update(app.config)

    # Ensure every task runs inside the Flask app-context
    class ContextTask(celery_app.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app.Task = ContextTask
    return celery_app


def create_app() -> Flask:
    app = Flask(__name__)

    # Flask (and Celery) configuration
    app.config.update(
        SECRET_KEY="your-secret-key",
        broker_url="redis://localhost:6379/0",
        result_backend="redis://localhost:6379/0",
    )

    # Attach Celery
    global celery  # bind to the module-level name
    celery = make_celery(app)

    return app
