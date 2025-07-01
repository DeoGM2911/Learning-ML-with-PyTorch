# app.py - main file for the polymer web app
#
# @author: Dung Tran
# @date: June 29, 2025
# @version: 1.0.0

from flask import Flask, render_template, request, jsonify
from app.tasks import add
from celery.result import AsyncResult

app = Flask(__name__)

@app.route('/', methods=['GET'])
def show_add_page():
    return render_template('add.html')

@app.route('/start-task', methods=['POST'])
def start_task():
    data = request.get_json()
    x = data.get('x')
    y = data.get('y')
    result = add.delay(x, y)
    return jsonify({"task_id": result.id})

@app.route('/status/<task_id>')
def get_status(task_id):
    result = AsyncResult(task_id)
    return jsonify({
        "status": result.status,
        "result": result.result if result.ready() else None
    })
