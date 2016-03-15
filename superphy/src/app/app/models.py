"""
models.py
"""
#pylint: disable=C0103, W0406

from passlib.apps import custom_app_context as pwd_context

from .. import db

class User(db.Model):
    """
    User model for database.
    """
    tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(32), index=True)
    password_hash = db.Column(db.String(128))

    def hash_password(self, password):
        """
        Takes a plain password as argument and stores a hash of it with the
        user
        """
        self.password_hash = pwd_context.encrypt(password)

    def verify_password(self, password):
        """
        ttakes a plain password as argument and returns True if the password
        is correct or False if not.
        """
        return pwd_context.verify(password, self.password_hash)

    def generate_auth_token(self, expiration=600):
        """
        the token is an encrypted version of a dictionary that has the id of
        the user. The token will also have an expiration time embedded in it,
        which by default will be of ten minutes (600 seconds).
        """
        s = Serializer(app.config['SECRET_KEY'], expires_in=expiration)
        return s.dumps({'id': self.id})

    @staticmethod
    def verify_auth_token(token):
        """
        If the token can be decoded then the id encoded in it is used to load
        the user, and that user is returned.
        """
        s = Serializer(app.config['SECRET_KEY'])
        try:
            data = s.loads(token)
        except SignatureExpired:
            return None # valid token, but expired
        except BadSignature:
            return None # invalid token
        user = User.query.get(data['id'])
        return user
