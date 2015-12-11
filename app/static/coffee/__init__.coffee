meta = new App()
home = new Home()

m.route(document.body, "/", {
    "/": home
    "/home": home
    "/meta": meta
    "/gbrowse": meta
    "/groups": meta
})