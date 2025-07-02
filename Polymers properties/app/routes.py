from flask import Blueprint, render_template, request, jsonify
from celery.result import AsyncResult
from . import celery
from .tasks import add

bp = Blueprint("main", __name__)

@bp.route("/")
def index():
    return render_template("index.html")

@bp.route("/add")
def show_add_page():
    return render_template("add.html")

@bp.route("/start-task", methods=["POST"])
def start_task():
    data = request.get_json()
    result = add.delay(data["x"], data["y"])
    return jsonify({"task_id": result.id})

@bp.route("/status/<task_id>")
def get_status(task_id):
    result = celery.AsyncResult(task_id)
    return jsonify({
        "status": result.status,
        "result": result.result if result.ready() else None
    })
