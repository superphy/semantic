var links;

links = [
  {
    title: 'item 1',
    url: '/item1'
  }
];

m.render(document.body, [
  m('ul.nav', [
    m('li', links.map(function(link) {
      return m('a', {
        href: link.url,
        config: m.route
      }, link.title);
    }))
  ])
]);


var getFoo;

var getEndpoint;

getEndpoint = function(url) {
  url = url || ""
  var response;
  console.log("http://" + location.hostname + "/" + url);
  response = m.request({
    method: "GET",
    url: "http://" + location.hostname + "/" + url,
    data: {},
    datatype: 'json',
    type: function(response) {
      localStorage.setItem(url, JSON.stringify(response));
      return response;
    }
  });
  return response;
};

m.render(document.body, [
	m('.', [
		getEndpoint('Foo')
	])
]);