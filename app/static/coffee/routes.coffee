class Routes
    locations: ->
        m.route(document.body, "/", {
          "/": home
          "/home": home
          #"/login": login
          "/hello": hello
      })

routes = new Routes()
routes.locations()





