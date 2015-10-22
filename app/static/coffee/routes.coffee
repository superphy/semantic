class Routes
    locations: ->
        m.route(document.body, "/", {
          "/": home
          "/home": home
          "/login": login
      })

routes = new Routes()
routes.locations()





