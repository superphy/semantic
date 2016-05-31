from flask import Flask, render_template, jsonify, json
import os
app = Flask(__name__)

@app.route("/Foo")
def hello():
	return jsonify({"HELLO":"WORLD"})

@app.route("/")
def foo():
	print os.path.abspath(os.path.dirname(__file__))
	return render_template('index.html')

@app.errorhandler(Exception)
def exception_handler(error):
    return jsonify({"ERROR": repr(error)})

if __name__ == "__main__":
    app.run()
