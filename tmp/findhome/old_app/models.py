from app import db

class User(db.Model):
	id = db.Column(db.Integer, primary_key=True)
	email = db.Column(db.String(254), index=True, unique=True)
	password = db.Column(db.String(128))
	first_name = db.Column(db.String(255))
	last_name = db.Column(db.String(255))
	org = db.Column(db.String(255))

	@property
	def is_authenticated(self):
		return True

	@property
	def is_active(self):
		return True

	@property #Nessesary?
	def is_anonymous(self):
		return False

	def get_id(self):
		return unicode(self.id)  # python 2
	def __repr__(self):
		return '<Email %r>' % (self.email)