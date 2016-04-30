class Page
    @view: (args) ->
        args = args or {}
        m('.', m.component(Navbar), args)