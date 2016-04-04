class Singleton
    model: () ->
    constructor: () ->
        @model()
    @get: ->
        @instance ?= new @()
        return @instance
    @getView: () ->
        return @get().view()
    view: () ->

class PageTemplate extends Singleton
    model: () ->
    controller: () ->
    view: () ->
        return [
            Header.getView()
        ]


class ComponentTemplate
    @controller: (args) ->
        view: () ->
            alert("Why are you calling a template view?")
    @view: (ctrl) ->
        return ctrl.view() #Should it pass in ctrl?

class Navbar extends ComponentTemplate
    @controller: (args) ->
        links = [
            {title: "Genome", url: "/Home.get()"}
            {title: "Group Browse", url: "/gbrowse"}
            {title: "Group Analyses", url: "/groups"}
            {title: "VF and AMR", url: "/factors"}
            {title: "Meta", url: "/meta"}
        ]
        view: () ->
            m("div", {class:'container-fluid'}, [
                m("nav", {class:'navbar navbar-default navbar-fixed-top',
                role:'navigation'}, [
                    m("div", {class:'navbar-header'}, [
                        m("a", {class:"", href:"/home", config:m.route},
                        "SuperPhy")
                    ])
                    m("ul", {class:'nav navbar-nav'}, [
                        m("li", m("a", {href: link.url, config:m.route},
                        link.title)) for link in links
                        m("li", {class:"dropdown"}, [
                            m("a", {
                                href:""
                                role:"button"
                                class:"dropdown-toggle"
                                'data-toggle':"dropdown"
                                }
                                "My Data", [
                                    m("b", {class:"caret"})
                                ]
                            )
                        ])
                    ])
                ])
            ])