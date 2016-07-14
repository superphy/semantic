
# coffeelint: disable=max_line_length

class Navbar
    @controller: (args) ->
        @links = [
            {title: "Genome", url: "/Home.get()"}
            {title: "Group Browse", url: "/gbrowse"}
            {title: "Group Analyses", url: "/groups"}
            {title: "VF and AMR", url: "/factors"}
            {title: "Meta", url: "/meta"}
            {title: "Test Table", url: "?/test"} # Shortcut to the table I'm working on
            {title: "New VF & AMR", url: "?/new"} # Not yet re-done
        ]
        return @
    @view: (ctrl) ->
        m(".container-fluid"
            m("nav.navbar navbar-default navbar-fixed-top"
                role:'navigation',
                m(".navbar-header"
                    m("a.navbar-brand", {href:"/home", config:m.route}, "SuperPhy")
                )
                m("ul.nav navbar-nav"
                    m("li"
                        m("a", {href: link.url, config:m.route}, link.title)
                    ) for link in ctrl.links
                    m("li.dropdown"
                        m("a.dropdown-toggle"
                            {href:"", role:"button", 'data-toggle':"dropdown"}
                            "My Data"
                            m("b.caret")
                        )
                    )
                )
            )
        )