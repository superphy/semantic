from flask import render_template, flash, redirect, session
from app import app
from .forms import LoginForm
from .forms import SignUpForm
from wtforms.validators import ValidationError

@app.route('/')
@app.route('/index')
def index():
    user = {'nickname': 'Miguel'}
    posts = [
        {
            'author': {'nickname': 'John'},
            'body': 'Beautiful day in Portland!'
        },
        {
            'author': {'nickname': 'Susan'},
            'body': 'The Avengers movie was so cool!'
        }
    ]
    return render_template('index.html',
                           title='Home',
                           user=user,
                           posts=posts)


@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        flash( "email: %s" % (form.email.data))
        flash( "password: %s" % (form.password.data))
        flash( "remember_me: %s" % (form.remember_me.data))
        return redirect('/index')
    return render_template('login.html',title='Log In',form=form,)

@app.route('/signup', methods=['GET', 'POST'])
def signup():
    form = SignUpForm()
    if form.validate_on_submit():
        flash( "email: %s" % (form.email.data))
        flash( "password: %s" % (form.password.data))
        flash( "first_name: %s" % (form.first_name.data))
        flash( "last_name: %s" % (form.last_name.data))
        flash( "org: %s" % (form.org.data))
        flash( "remember_me: %s" % (form.remember_me.data))
        flash( "You have successfully signed up! An email has been sent to %s." % (form.email.data))
        return redirect('/index')
    return render_template('signup.html',title='Sign Up',form=form,)

@app.errorhandler(404)
def page_not_found(e):
	return render_template('404.html'), 404

@app.errorhandler(500)
def internal_server_error(e):
    return render_template('500.html'), 500

'''@lm.user_loader
def load_user(id):
    return User.query.get(int(id))'''
