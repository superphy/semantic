class GeneForm
    constructor: (@name, @type) ->
        @model()
        @controller()


    model: =>
        @data ?= new GeneData(@type)
        @genelist ?= new GeneList()

        @searchterm = m.prop('')


    controller: =>
        escapeRegExp = (str) -> str.replace(/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]/g, "\\$&")

        @search = (term) =>
            @searchterm = m.prop(term)
            regex = new RegExp(escapeRegExp(term), "i")

            for row in @data.rows
                gene = row.Gene_Name

                if regex.test(gene)
                    row.visible(true)
                else
                    row.visible(false)

            return true

        @select = (all) =>
            for row in @data.rows
                if all is "true" then row.selected(true) else row.selected(false)

        @toggle = (checked) =>
            console.log(checked)
            if checked is "true" then row.selected(false) else row.selected(true)
    

    view: () =>
        gene_name = @name.toLowerCase()

        m('div', {class: 'panel panel-default'}, [
            m('div', {class: 'panel-heading', id: 'vf-panel-header'}, [
                m('h4', {class: 'panel-title'}, m('a', {href: '#vf-form', config:m.route}, "#{@name} Form"))
            ])
            m('div', {class: 'panel', id: 'vf-panel'}, [
                m('div', {class: 'panel-body'}, [
                    m('div', {class: 'row'}, [
                        m('div', {class: 'col-md-6 col-md-offset-3'}, [
                            m('div', {class: 'selected-gene-list-wrapper', id: 'vf-selected-list'}, [
                                m('fieldset', [
                                    m('span', ['Selected factors:'])
                                    m('ul', {id: 'vf-selected'}, [
                                        for row in @data.rows
                                            if row.selected()
                                                m('li', {class: 'selected-gene-item'}, [
                                                    row["Gene_Name"]
                                                    # m('a', {href: "#", checked: row.selected(), onclick: m.withAttr("checked", toggle)}, [
                                                    #     m('i', {class: 'fa fa-times'})
                                                    # ])
                                                ])
                                    ])
                                ])
                            ])
                        ])
                    ])

                    m('div', {class: 'row'}, [
                        m('div', {class: 'gene-search-control-row'}, [
                            m('div', {class: 'col-md-3'}, [
                                m('input[type=text]', {id: 'vf-autocomplete', \
                                                       class: 'form-control', \
                                                       placeholder: "Search #{gene_name} gene in list", \
                                                       value: @searchterm(), \
                                                       onkeyup: m.withAttr("value", @search)})
                            ])
                            m('div', {class: 'col-md-3'}, [
                                m('div', {class: 'btn-group'}, [
                                    m('button', {class: 'btn btn-link', checked: true, \
                                                 onclick: m.withAttr("checked", @select)}, "Select All")
                                    m('button', {class: 'btn btn-link', checked: false, \
                                                 onclick: m.withAttr("checked", @select)}, "Deselect All")
                                ])
                            ])
                        ])
                    ])

                    m('div', {class: 'row'}, [
                        m('div', {class: 'cold-md-6'}, [
                            m('div', {class: 'gene-list-wrapper'}, [
                                m('fieldset', [
                                    m('span', {class: 'col-md-12'}, ["Select one or more #{@type} factors"])
                                    m('div', {class: 'col-md-12'}, [
                                        m('div', {class: 'superphy-table', id: 'vf-table'}, [
                                            @genelist.view(@data)
                                        ])
                                    ])
                                ])
                            ])
                        ])
                    ])

                    m('div', {class: 'row'}, [
                        m('div', {class: 'cold-md-6'}, [
                            m('div', {class: 'gene-category-wrapper'}, [
                                m('div', {class: 'gene-category-intro'}, [
                                    m('span', "Select category to refine list of genes:")
                                    for category in @data.categories
                                        m('div', {class: "row"}, [
                                            m('div', {class: "col-xs-12"}, [
                                                m('hello')
                                            ])
                                        ])
                                ])
                            ])
                        ])
                    ])
                ])
            ])
        ])