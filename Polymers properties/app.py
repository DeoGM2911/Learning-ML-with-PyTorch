# app.py - main file for the polymer web app
#
# @author: Dung Tran
# @date: June 29, 2025
# @version: 1.0.0

from flask import Flask

app = Flask(__name__)

@app.route('/')
def hello_world():
    return 'Hello, World!'

if __name__ == '__main__':
    app.run(debug=True)