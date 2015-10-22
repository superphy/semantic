authenticate ={
  controller: (p1, p2) ->
      runAuth: (p1, p2) ->

        m.request({
          method:"POST"
          url:"scripts/scrypt_pass.py?params=" + p1 + "," + p2
        })

  view: (ctrl) ->
    m('div', ctrl.runAuth())
}


#var returnValue = m.prop([]);
#
#db.insert = function(tableName){
#  this.rValue = m.request({method: "POST"
#    ,url:"scripts/db.pl?mode=create&params=" + tableName
#  });
#};
#
#db.controller = function(){
#  this.tableName = m.route.param("tableName");
#db.insert(this.tableName);
#};
#
#db.view = function(){
#return m("div","DONE");
#}
