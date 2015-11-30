class Routes
      m.route(document.body, "/", {
          "/": home
          "/home": home
          "/meta": new Page(limit = 20, page = 2)
      })

