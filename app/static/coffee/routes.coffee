class Routes
      m.route(document.body, "/", {
          "/": home
          "/home": home
          "/meta": Meta_Data.get()
      })

