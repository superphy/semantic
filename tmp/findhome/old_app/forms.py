from flask.ext.wtf import Form
from wtforms import StringField, BooleanField, PasswordField
from wtforms.validators import DataRequired
from wtforms.validators import Email
from wtforms.validators import Optional
from wtforms.validators import ValidationError

_debugging = True

if _debugging:
    email_ = "BryceDrew.1@gmail.com"
    first_name_ = "Bryce"
    last_name_ = "Drew"
    org_ = "Laboratory for Foodborne Zoonoses Canada"
    remember_me_ = True
else:
    email_ = ""
    first_name_ = ""
    last_name_ = ""
    org_ = ""
    remember_me_ = False

def unique_email(form, email):
    return True
    #except Email.AlreadyExists:
    #   raise ValidationError("%s is already registered" % email)

#Include all validation of fields as methods in the validators section.
class LoginForm(Form):
    email = StringField('Email', validators=[DataRequired(),Email()], default=email_)
    password = PasswordField('Password', validators=[DataRequired()])
    remember_me = BooleanField('remember_me', default=remember_me_)

class SignUpForm(Form):

    email = StringField('Email', validators=[DataRequired(),Email(),unique_email], default=email_)
    password = PasswordField('Password', validators=[DataRequired()])
    first_name = StringField('First_name',validators=[Optional()], default=first_name_)
    last_name = StringField('Last_name',validators=[Optional()], default=last_name_)
    org = StringField('Organization',validators=[Optional()], default=org_)
    remember_me = BooleanField('remember_me', default=remember_me_)
    




#Include Genome submission and querying forms here