###
CLASS GeneSearchPanel

Component class for the Gene Search Panel for one gene type
  in the VF and AMR gene feature.

Args passed in:
    title: title of the panel
    type: string that describes type of gene (either "vf" or "amr")
    data: data model (global model that holds the selected genes & genomes)
    categories: object of categories for the genes
###

class GeneSearchPanel
    @controller: (args) ->
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
                            data: ctrl.data
                            categories: ctrl.categories
                        })
                    ])
                ])
            ])
        ]) 


###
COMPONENT SelectedGenes

Component that displays the selected genes in the search search panel.

Args passed in:
    data: data model (global model that holds the selected genes & genomes)
    type: string that describes type of gene (either "vf" or "amr")
###
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


###
COMPONENT SearchSelect

Search bar for the gene table.

Args passed in:
    data: data model (global model that holds the selected genes & genomes)
    type: string that describes type of gene (either "vf" or "amr")
###
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


###
COMPONENT GeneTable

Component that displays the gene table.

Args passed in:
    data: data model (global model that holds the selected genes & genomes)
    type: string that describes type of gene (either "vf" or "amr")
###
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


###
COMPONENT GeneTable

Component for the category selection (whole section, not just one select bar)

Args passed in:
    data: data model (global model that holds the selected genes & genomes)
    categories: object that has the categories and subcategories for the
                particular set of genomes.
###
CategorySelection =
    controller: (args) ->
        @count = 0
        @data = [
            {id: 1, name: "John"},
            {id: 2, name: "Mary"},
            {id: 3, name: "Seniqua"}
          ]
        @currentUser = m.prop("Select Category")
        @changeUser = (id) ->
            console.log(id)
        return @
    view: (ctrl, args) ->
        m(".col-md-6", [
            m(".gene-category-wrapper", [
                m(".gene-category-intro", [
                    m('span', "Select category to refine list of genes:")
                ])
                for category of args.categories
                    [m('div', {class: "row"}, [
                        m('div', {class: "category-header col-xs-12"}, [
                            category
                        ])
                    ])
                    # m('div', {class: "selectize-control form-control single"}, [
                    #     m('div', {class: "selectize-input items not-full has-options"}, [
                    #         m('input', {type: "text", autocomplete: "off", placeholder: "--Select a category--", style: "width: 134px"})
                    #     ])
                    # ])

                    # m('select', {class: "subcategories", multiple: "multiple"}, [
                    #     for subcat in args.categories[category]
                    #         m('option', subcat)
                    # ])
                    m('.row',
                        m('.col-xs-12',
                            # m('.', {class: "selectize-control form-control single"}, [
                            #     m('div', {class: "selectize-input items not-full has-options"}, [
                            m.component(SelectMult, {
                                #filterbar: new SelectMult()
                                category: category
                                data: args.categories[category]
                                value: ctrl.currentUser
                                onchange: ctrl.changeUser
                            })
                            #     ])
                            # ])
                        )
                    )
                    ]
            ])
        ])

SelectBar = 
    controller: (args) ->
        @filterCategories = (category, subcategory) ->
            console.log("Filtering for ", category, " genes")
            for row in args.data.rows
                cat = row.Category
                subcat = row.Sub_Category
                if cat is category and subcat is subcategory
                    row.visible(true)
                else
                    row.visible(false)

        @onchange = (cat, sub) ->
            @filterCategories(cat, sub)

        return @

    config: (ctrl) ->
        return (element, isInitialized) ->
            el = $(element)
            if not isInitialized
                el.on("change", (e) ->
                    if typeof ctrl.onchange is "function"
                        ctrl.onchange()
                )
    view: (ctrl, args) ->
        m('select', {"data-category-id": args.category, multiple: "multiple", class:"form-control"}, [
            for subcat,index in args.data
                attrs = {id: "#{index}", value: subcat, title:subcat}
                m('option', attrs, subcat)

            m('option', {value:'null', selected: "selected"}, "--Select Category--")
        ])

SelectMult =
    config: (ctrl) ->
        return (element, isInitialized) ->
            el = $(element)
            if not isInitialized
                attrs = {
                    placeholder: "--Select a category--"
                }
                el.select2(attrs)
            #         .on("change", (e) ->
            #             m.startComputation()
            #             if typeof ctrl.onchange is "function"
            #                 ctrl.onchange(el.select2("val"))
            #                 #ctrl.filter(args.category, item)

            #             m.endComputation()
            #         )
    view: (ctrl, args) ->
        #selectedId = args.value()
        return m("select", {"data-category-id": args.category, \
                            multiple: "multiple", \
                            class: "form-control", \
                            config: SelectMult.config(args)}, [
            for subcat,index in args.data
                attrs = {id: "#{index}", value: subcat, title:subcat}
                m('option', attrs, subcat)

            #m('option', {value:'null', selected: "selected"}, "--Select Category--")
        ])

Select2 =
    controller: (args) ->
        @filter = (category, subcategory) ->
            console.log("Filtering for ", category, " genes")
            for row in args.data.rows
                cat = row.Category
                subcat = row.Sub_Category
                if cat is category and subcat is subcategory
                    row.visible(true)
                else
                    row.visible(false)

        return @

    config: (ctrl) ->
        return (element, isInitialized) ->
            if typeof jQuery isnt 'undefined' and typeof jQuery.fn.select2 isnt 'undefined'
                el = $(element)
                if not isInitialized
                    el.select2()
                        .on("change", (e) ->
                            #id = el.select2("val")
                            m.startComputation()
                            # ctrl.data.map((d) ->
                            #     if d is id
                            #         ctrl.value(d)
                            # )

                            if typeof ctrl.onchange is "function"
                                ctrl.onchange(el.select2("val"))
                                #ctrl.filter(args.category, item)

                            m.endComputation()
                        )

                el.select2("val", ctrl.value)
                #el.val(ctrl.value()).trigger("change")
            else
                console.warn('ERROR: You need jquery and Select2 in the page')

    view: (ctrl, args) ->
        #selectedId = args.value()
        return m("select", {config: Select2.config(args)}, [
            args.data.map((item) ->
                attrs = {value:item}
                # if item is selectedId
                #     attrs.selected = "selected"
                return m("option", attrs, item) )
        ])