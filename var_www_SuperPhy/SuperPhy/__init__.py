
import os
from flask import Flask, render_template, jsonify

from SuperPhy.blueprints.simple import simple as simple_blueprint

app = Flask(__name__)

@app.route("/")
def index():
	return render_template('simple_index.html')

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

app.register_blueprint(simple_blueprint, url_prefix='/simple')


if __name__ == "__main__":
    app.run()

"""
    from.example import example as example_blueprint
    app.register_blueprint(example_blueprint, url_prefix='/example')

    from .data import data as data_blueprint
    app.register_blueprint(data_blueprint, url_prefix='/data')

    from .api import api as api_blueprint
    app.register_blueprint(api_blueprint, url_prefix='/api')

    from .upload import upload as upload_blueprint
    app.register_blueprint(upload_blueprint, url_prefix='/upload')

    return app
"""