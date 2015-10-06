from flask.ext.wtf import Form
from wtforms import StringField, BooleanField, PasswordField
from wtforms.validators import DataRequired
from wtforms.validators	import Email
from wtforms.validators	import Optional
from wtforms.validators	import ValidationError

#move to validation
def unique_email(form, email):
	#Check if email is unique
	raise ValidationError("This email is already registered")

#Include all validation of fields as methods in the validators section.
class LoginForm(Form):
    email = StringField('Email', validators=[DataRequired(),Email()])
    password = PasswordField('Password', validators=[DataRequired()])
    remember_me = BooleanField('remember_me', default=False)

class SignUpForm(Form):
    email = StringField('Email', validators=[DataRequired(),Email(),unique_email])
    password = PasswordField('Password', validators=[DataRequired()])
    first_name = StringField('First_name',validators=[Optional()])
    last_name = StringField('Last_name',validators=[Optional()])
    org = StringField('Organization',validators=[Optional()])
    remember_me = BooleanField('remember_me', default=False)






#Include Genome submission and querying forms here