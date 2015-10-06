CSRF_ENABLED = True
WTF_CSRF_ENABLED = True
SECRET_KEY = 'you-will-never-guess'
#https://pythonhosted.org/Flask-Security/configuration.html
# Specifies the name for the Flask-Security blueprint. Defaults to security.
#SECURITY_BLUEPRINT_NAME
# Specifies the URL prefix for the Flask-Security blueprint. Defaults to None.			
#SECURITY_URL_PREFIX
# Specifies whether or not to flash messages during security procedures. Defaults to True.
#SECURITY_FLASH_MESSAGES
# Specifies the password hash algorithm to use when encrypting and decrypting passwords. Recommended values for production systems are bcrypt, sha512_crypt, or pbkdf2_sha512. Defaults to plaintext.
SECURITY_PASSWORD_HASH = 'bcrypt'
# Specifies the HMAC salt. This is only used if the password hash type is set to something other than plain text. Defaults to None.
SECURITY_PASSWORD_SALT = SECRET_KEY
# Specifies the email address to send emails as. Defaults to no-reply@localhost.			
#SECURITY_EMAIL_SENDER
# Specifies the query sting parameter to read when using token authentication. Defaults to auth_token.		
#SECURITY_TOKEN_AUTHENTICATION_KEY
# Specifies the HTTP header to read when using token authentication. Defaults to Authentication-Token.
#SECURITY_TOKEN_AUTHENTICATION_HEADER
# Specifies the default authentication realm when using basic HTTP auth. Defaults to Login Required
#SECURITY_DEFAULT_HTTP_AUTH_REALM

#https://pythonhosted.org/Flask-Security/models.html 
#If you enable account confirmation by setting your application's 
#SECURITY_CONFIRMABLE configuration value to True, your User model 
#will require the following additional field: confirmed_at
#SECURITY_CONFIRMABLE = True

#https://pythonhosted.org/Flask-Security/models.html
#If you enable user tracking by setting your application's 
#SECURITY_TRACKABLE configuration value to True, your User 
#model will require the following additional fields: 
#    last_login_at
#    current_login_at
#    last_login_ip
#    current_login_ip
#    login_count'''
SECURITY_TRACKABLE = True