using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;

namespace EarlyBirdnBird.Pages.Host
{
    public class ListingsModel : PageModel
    {

        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;
        public IList<Rental> Rentals { get; set; }


        public ListingsModel(ApplicationDbContext db, UserManager<AppUser> um)
        {
            _dbContext = db;
            _userManager = um;
        }

        public async Task<IActionResult> OnGetAsync()
        {
            var currentUser = await _userManager.GetUserAsync(User);

            if (currentUser == null)
            {
                // User not found, handle the error accordingly
                return NotFound("User not found.");
            }

            var rentals = from r in _dbContext.Rentals 
                          where r.OwnerName == currentUser
                          select r;


            Rentals = await rentals
                .Include(r => r.Images)
                .Include(r => r.Prices)
                .Include(r => r.OwnerName)
            .ToListAsync();
            return Page();
        }
    }
}
