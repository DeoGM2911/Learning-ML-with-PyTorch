from flask import Blueprint, render_template, request, jsonify
from .tasks import process_prediction

bp = Blueprint("main", __name__)

@bp.route("/")
def index():
    return render_template("index.html")


@bp.route("/predictions")
def predictions():
    return render_template("predictions.html")


@bp.route('/api/predict', methods=['POST'])
def predict():
    text_input = request.form.get('text')
    file_obj   = request.files.get('file')

    if not text_input and not file_obj:
        return jsonify({'error': 'No input provided.'}), 400
    
    # read file into bytes if present
    file_bytes = file_obj.read() if file_obj else None
    
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
            'reply_html': task.result
        })
    else:
        # FAILURE, RETRY, etc
        return jsonify({
            'status': task.state,
            'error': task.result
        }), 500
