class GeneSearchPanel
    @controller: (args) ->
        @toggle = (checked) ->

        return @

    @view: (args, ctrl) ->
        m(".panel panel-default", [
            m(".panel-heading", {id: "#{ctrl.type}-panel-header"}, [
                m("h4", {class: "panel-title"},
                    m("a", {href: "#vf-form", config:m.route}, "#{ctrl.title} Form"))
            ])
            m(".panel", {id: "#{ctrl.type}-panel"}, [
                m(".panel-body", [
                    m(".row", [
                        m.component(SelectedGenes, {
                            data: ctrl.data
                            type: ctrl.type
                        })
                    ])
                    m(".row", [
                        m.component(SearchSelect, {
                            title: ctrl.title
                            type: ctrl.type
                            data: ctrl.data
                        })
                    ])
                    m(".row", [
                        m.component(GeneTable, {
                            data: ctrl.data
                            type: ctrl.type
                        })
                        m.component(CategorySelection, {
                            ## Test values
                            categories: []
                            
                        })
                    ])
                ])
            ])
        ]) 


SelectedGenes =
    controller: (args) ->
    view: (ctrl, args) ->
        m(".", {class: "col-md-6 col-md-offset-3"}, [
            m(".", {class: "selected-gene-list-wrapper", \
                    id: "#{args.type}-selected-list"}, [
                m("fieldset", [
                    m("span", ["Selected factors:"])
                    m("ul", {id: "#{args.type}-selected"}, [
                        for row in args.data.rows when row.selected()
                            m('li', {class: 'selected-gene-item'}, [
                                row["Gene_Name"]
                            ])
                    ])
                ])
            ])
        ])

SearchSelect =
    controller: (args) ->
        self = @
        @searchterm = m.prop("")
        escapeRegExp = (str) ->
            str.replace(/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]/g, "\\$&")

        @search = (term) ->
            self.searchterm = m.prop(term)
            regex = new RegExp(escapeRegExp(term), "i")

            console.log("Searching.....")

            for row in args.data.rows
                gene = row.Gene_Name
                if regex.test(gene)
                    row.visible(true)
                else
                    row.visible(false)
            return true

        @select = (all) ->
            for row in args.data.rows
                if all is "true" then row.selected(true) else row.selected(false)

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
                                 onclick: m.withAttr("checked", ctrl.select)}, \
                                 "Select All")
                    m('button', {class: 'btn btn-link', checked: false, \
                                 onclick: m.withAttr("checked", ctrl.select)}, \
                                 "Deselect All")
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
                        m(".",
                            m.component(TableView, {
                                data: args.data
                                checkbox: true
                            })
                        )
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
