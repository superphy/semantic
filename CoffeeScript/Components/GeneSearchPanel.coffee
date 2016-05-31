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
# coffeelint: disable=max_line_length

class GeneSearchPanel
    ###
    COMPONENT SelectedGenes

    Component that displays the selected genes in the search search panel.

    Args passed in:
        data: data model (global model that holds the selected genes & genomes)
        type: string that describes type of gene (either "vf" or "amr")
    ###
    SelectedGenes = (args) ->
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
    COMPONENT GeneTable

    Component that displays the gene table.

    Args passed in:
        data: data model (global model that holds the selected genes & genomes)
        type: string that describes type of gene (either "vf" or "amr")
    ###
    
    GeneTable = (args) ->
        m(".col-md-6", [
            m(".gene-list-wrapper", [
                m("fieldset", [
                    m("span.col-md-12", [
                        "Select one or more #{args.type} factors"
                        m(".",
                            m.component(Table,
                                data: args.data
                                checkbox: true
                            )
                        )
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
    class SearchSelect
        @controller: (args) ->
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

            return @

        @view: (ctrl, args) ->
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
    COMPONENT CategorySelection

    Component for the category selection (whole section, not just one select bar)

    Args passed in:
        data: data model (global model that holds the selected genes & genomes)
        categories: object that has the categories and subcategories for the
                    particular set of genomes.
    ###
    class CategorySelection
        @controller: (args) ->
            @clear= () ->
                ## Clears filter of categories
                ## Bug: categories don't disappear in select boxes yet - look into
                ##     select2 documentation for clear under programmatic access :)
                for row in args.data.rows
                    row.visible(true)
            return @
        @view: (ctrl, args) ->
            m(".col-md-6", [
                m(".gene-category-wrapper", [
                    m(".gene-category-intro", [
                        m('span', "Select category to refine list of genes:")
                        m('button', {onclick: ctrl.clear}, "Clear")
                    ])
                    for category of args.categories
                        [m('.row', [
                            m('.category-header col-xs-12', [
                                category
                            ])
                        ])
                        m('.row',
                            m('.col-xs-12',
                                m.component(SelectBar, {
                                    categories: args.categories
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

    WORK TO DO:
    - Fix problem with selecting multiple genes from different main categories
      Currently, if I pick subcategories from say "Adherence", then pick
      some others from "Autotransporter", the "Adherence" ones will disappear.
    ###
    class SelectBar
        @controller: (args) ->
            @allcategories = {}
            @createCategories = () ->
                #console.log(args.categories)
                listCats = Object.keys(args.categories)
                for cat in listCats
                    ## boolean for indicating if the category has been selected
                    self.allcategories[cat] = false

            @filterCategories = (category, subcategories) ->
                console.log("Filtering for ", category, " genes")
                for row in args.data.rows
                    cat = row.Category
                    subcat = row.Sub_Category
                    if !subcategories or (cat is category and subcat in subcategories)
                        row.visible(true)
                    else
                        row.visible(false)
            return @

        @config: (ctrl, args) ->
            return (element, isInitialized) ->
                ctrl.createCategories()
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
        @view: (ctrl, args) ->
            return m("select", {"data-category-id": args.category, \
                                multiple: "multiple", \
                                class: "form-control", \
                                config: SelectBar.config(ctrl, args)}, [
                for subcat, index in args.subcategories
                    attrs = {id: "#{index}", value: subcat, title:subcat}
                    m('option', attrs, subcat)
            ])

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
                        SelectedGenes(
                            data: ctrl.data
                            type: ctrl.type
                        )
                    ])
                    m(".row", [
                        m.component(SearchSelect, {
                            type: ctrl.type
                            data: ctrl.data
                        })
                    ])
                    m(".row", [
                        GeneTable(
                            data: m.prop(ctrl.data)
                            type: ctrl.type
                        )
                        m.component(CategorySelection, {
                            data: ctrl.data
                            categories: ctrl.categories
                        })
                    ])
                ])
            ])
        ])
