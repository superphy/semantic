from superphy.sparql import user as sparql
from superphy import endpoint

'''
from superphy.user import *
email = "dummy"
insert_user(email)
'''

def create_user(email):
	#See if email exists
	if endpoint.ask(sparql.email_exists(email)):
		#This doesn't acutally stop anything with the sql though, so don't think this is a good place to put a validator. -Bryce
		print "Already Exists!"
		return False
	#Get next user id
	results = endpoint.query(sparql.last_user())
	user_id = 0
	for result in results["results"]["bindings"]:
		for i, thing in enumerate('s'):
			user_id = int(result[thing]['value'][28:-5]) + 1
			break
		break
	#email = "Test.%s@test.com" % user_id #This is here to provide an easy sample email when debugging.
	user_id = str(user_id).zfill(10)
	#Insert into graph
	endpoint.update(sparql.insert_user(user_id, email))