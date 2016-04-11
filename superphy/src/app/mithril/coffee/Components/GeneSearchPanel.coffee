class GeneSearchPanel
    @controller: (args) ->
        self = @
        @args = args || {}

        @select = (all) ->

        @toggle = (checked) ->

        return @

    @view: (ctrl) ->
        m(".panel panel-default", [
            m(".panel-heading", {id: "#{ctrl.args.type}-panel-header"}, [
                m("h4", {class: "panel-title"},
                    m("a", {href: "#vf-form", config:m.route}, "#{ctrl.args.title} Form"))
            ])
            m(".panel", {id: "#{ctrl.type}-panel"}, [
                m(".panel-body", [
                    m(".row", [
                        m.component(SelectedGenes, {
                            data: "data goes here eventually"
                            type: ctrl.args.type
                        })
                    ])
                    m(".row", [
                        m.component(SearchSelect, {
                            title: ctrl.args.title
                            type: ctrl.args.type
                        })
                    ])
                    m(".row", [
                        m.component(GeneTable, {
                            type: ctrl.args.type
                        })
                        m.component(CategorySelection)
                    ])
                ])
            ])
        ])


SelectedGenes =
    view: (ctrl, args) ->
        m(".", {class: "col-md-6 col-md-offset-3"}, [
            m(".", {class: "selected-gene-list-wrapper", id: "#{args.type}-selected-list"}, [
                m("fieldset", [
                    m("span", ["Selected factors:"])
                    m("ul", {id: "#{args.type}-selected"}, [
                        args.data
                    ])
                ])
            ])
        ])

SearchSelect =
    controller: (args) ->
        @searchterm = m.prop("")
        escapeRegExp = (str) ->
            str.replace(/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]/g, "\\$&")

        @search = (term) ->
            self.searchterm = m.prop(term)
            regex = new RegExp(escapeRegExp(term), "i")

            console.log("Searching.....")

            ## stuff about data
        return

    view: (ctrl, args) ->
        m(".gene-search-control-row", [
            m(".col-md-3", [
                m("input[type=text]", {id: "#{args.type}-autocomplete", \
                                       class: "form-control", \
                                       placeholder: "Search gene in list", \
                                       value: ctrl.searchterm(), \
                                       onkeyup: m.withAttr("value", ctrl.search)})
            ])
            m('.col-md-3', [
                m('.btn-group', [
                    m('button', {class: 'btn btn-link', checked: true, \
                                 onclick: m.withAttr("checked", @select)}, "Select All")
                    m('button', {class: 'btn btn-link', checked: false, \
                                 onclick: m.withAttr("checked", @select)}, "Deselect All")
                ])
            ])
        ])

GeneTable =
    controller: (args) ->
    view: (ctrl, args) ->
        m(".col-md-6", [
            m(".gene-list-wrapper", [
                m("fieldset", [
                    m("span", {class: "col-md-12"}, [
                        "Select one or more #{args.type} factors"
                    ])
                ])
            ])
        ])


CategorySelection =
    controller: (args) ->
    view: (ctrl, args) ->
        m(".col-md-6", [
            m(".gene-category-wrapper", [
                m(".gene-category-intro", [
                    m('span', "Select category to refine list of genes:")
                ])
            ])
        ])
