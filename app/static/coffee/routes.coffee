header = new Header()

class Routes
	home = new Home
	m.route(document.body, "/", {
		"/": home
		"/home": home
		"/meta": new Page(limit = 20, page = 2)
		"/gbrowse": new App
	})