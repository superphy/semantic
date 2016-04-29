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
        m(".col-md-6 col-md-offset-3", [
            m(".selected-gene-list-wrapper", \
                    id: "#{args.type}-selected-list", [
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
                    m("span.col-md-12", [
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
COMPONENT CategorySelection

Component for the category selection (whole section, not just one select bar)

Args passed in:
    data: data model (global model that holds the selected genes & genomes)
    categories: object that has the categories and subcategories for the
                particular set of genomes.
###
CategorySelection =
    controller: (args) ->
        return @
    view: (ctrl, args) ->
        m(".col-md-6", [
            m(".gene-category-wrapper", [
                m(".gene-category-intro", [
                    m('span', "Select category to refine list of genes:")
                ])
                for category of args.categories
                    [m('.row', [
                        m('.category-header col-xs-12', [
                            category
                        ])
                    ])
                    m('.row',
                        m('.col-xs-12',
                            m.component(SelectMult, {
                                category: category
                                subcategories: args.categories[category]
                                data: args.data
                            })
                        )
                    )
                    ]
            ])
        ])


###
COMPONENT SelectMult

Component for each category search box. Based on the Select2 jQuery library

Args:
    category: top category
    subcategories: list of subcategories
    data: gene search model
###
SelectMult =
    all_subcategories: []
    controller: (args) ->
        @filterCategories = (category, subcategories) ->
            console.log("Filtering for ", category, " genes")
            console.log("selectedcats:", subcategories)
            # for sc in SelectMult.all_subcategories
            #     console.log("the subcat is:", sc)
            for row in args.data.rows
                cat = row.Category
                subcat = row.Sub_Category
                if !subcategories or (cat is category and subcat in subcategories)
                    row.visible(true)
                else
                    row.visible(false)

        return @

    config: (ctrl, args) ->
        return (element, isInitialized) ->
            el = $(element)
            if not isInitialized
                attrs = {
                    placeholder: "--Select a category--"
                }
                el.select2(attrs)
                    .on("change", (e) ->
                        m.startComputation()
                        selected = []
                        if el.val() then selected = el.val()
                        ctrl.filterCategories(args.category, el.val())

                        m.endComputation()
                    )
    view: (ctrl, args) ->
        return m("select", {"data-category-id": args.category, \
                            multiple: "multiple", \
                            class: "form-control", \
                            config: SelectMult.config(ctrl, args)}, [
            for subcat, index in args.subcategories
                attrs = {id: "#{index}", value: subcat, title:subcat}
                m('option', attrs, subcat)
        ])