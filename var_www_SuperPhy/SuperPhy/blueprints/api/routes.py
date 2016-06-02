from flask import render_template

from SuperPhy import app
from SuperPhy.blueprints.api import api

@api.route('/')
def show_routes():
    """
    This function displays a template that shows to the user what Flask has
    registered as urls. This is for debugging, and illustration purposes.
    """
    routes = []
    for rule in app.url_map.iter_rules():
        if "GET" in rule.methods:
            routes.append(rule.rule)
    routes.sort()
    return render_template('routes.html', routes=sorted(routes, key=lambda s: s.lower()))
