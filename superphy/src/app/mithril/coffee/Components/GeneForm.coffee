class GeneForm
    constructor: (@name, @type) ->
        @model()

    model: =>
        @data ?= new GeneData(@type)
        @genelist ?= new GeneList()

        ## Array for selected genes for this form
        @selected = []

    controller: =>

    
    view: () =>
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
                                            if row["selected"]()
                                                m('li', {class: 'selected-gene-item'}, [
                                                    row["Gene_Name"]])
                                    ])
                                ])
                            ])
                        ])
                    ])

                    m('div', {class: 'row'}, [
                        m('div', {class: 'gene-search-control-row'}, [
                            m('div', {class: 'col-md-3'}, [
                                m('input', {class: 'form-control', id:'vf-autocomplete', placeholder: "Search #{@name} gene in list"})
                            ])
                            m('div', {class: 'col-md-3'}, [
                                m('input', {class: 'btn-group'}, [
                                    m('button', {class: 'btn btn-link', id: 'vf-select-all'}, ["Select All"])
                                    m('button', {class: 'btn btn-link', id: 'vf-deselect-all'}, ["Deselect All"])
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
                ])
            ])
        ])