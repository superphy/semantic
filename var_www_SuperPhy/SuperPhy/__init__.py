
import os
from flask import Flask, render_template, jsonify

app = Flask(__name__)

@app.route("/")
def index():
	print os.path.abspath(os.path.dirname(__file__))
	return render_template('index.html')

@app.route("/Info")
def routes():
    routes = []
    for rule in app.url_map.iter_rules():
        if "GET" in rule.methods:
            routes.append(rule.rule)
    routes.sort()
    return render_template('routes.html', routes=sorted(routes, key=lambda s: s.lower()))

@app.errorhandler(Exception)
def exception_handler(error):
    return jsonify({"ERROR": repr(error)})

if __name__ == "__main__":
    app.run()
