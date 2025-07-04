from flask import Blueprint, render_template, request, jsonify
from celery.result import AsyncResult
from . import celery
from .tasks import add, process_prediction
from PIL import Image
from services.images_to_smiles.utils.validate_txt import clean_smiles, is_smiles, is_smiles_polymer

bp = Blueprint("main", __name__)

@bp.route("/")
def index():
    return render_template("index.html")


@bp.route('/api/predict', methods=['POST'])
def predict():
    text_input = request.form.get('text')
    file_obj   = request.files.get('file')

    if not text_input and not file_obj:
        return jsonify({'error': 'No input provided.'}), 400
    
    # read file into bytes if present
    file_bytes = file_obj.read() if file_obj else None

    # validate the text input
    if not is_smiles(text_input):
        return jsonify({'error': 'Invalid SMILES input.'}), 400
    text_input = clean_smiles(text_input)
    if not is_smiles_polymer(text_input):
        return jsonify({'error': 'Input is not a valid polymer.'}), 400
    
    # enqueue the Celery task
    task = process_prediction.delay(text_input, file_bytes)

    # immediately return the task ID
    return jsonify({
        'task_id': task.id,
        'status':  'pending'
    }), 202


@bp.route('/api/result/<task_id>', methods=['GET'])
def get_result(task_id):
    task = process_prediction.AsyncResult(task_id)
    if task.state == 'PENDING':
        return jsonify({'status': 'pending'}), 202
    elif task.state == 'SUCCESS':
        return jsonify({
            'status': 'done',
            'reply':  task.result
        })
    else:
        # FAILURE, RETRY, etc
        return jsonify({
            'status': task.state,
            'error': task.result
        }), 500



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
