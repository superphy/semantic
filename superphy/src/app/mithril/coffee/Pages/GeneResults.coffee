# coffeelint: disable=max_line_length

class GeneResults extends Page
    @model: () ->
        @title = "Virulence Factor and AMR Results"
        return @

    @controller: (args) ->
        @title = GeneResults.model()
        return @
    @view: (ctrl) ->
        super m('.'
            m('div', {class: 'toc'}, [
                m('div', {class: 'well toc-well'}, [
                    m('div', {id: 'vf_result_legend', class: 'legend', display: 'none'}, [
                        ctrl.title
                    ])
                ])
            ])
            m('div', {id: 'results'}, [
                m('div', {id: 'vf_results'},[
                    m('hr')
                    m('h4', "Detected Virulence Factor Alleles")
                ])
            ])
        )

