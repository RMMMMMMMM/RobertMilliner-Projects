using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;

namespace EarlyBirdnBird.Pages.Host
{
    public class ReviewsModel : PageModel
    {
        private readonly UserManager<AppUser> _userManager;
        private readonly ApplicationDbContext _dbContext;

        public ReviewsModel(UserManager<AppUser> userManager, ApplicationDbContext dbContext)
        {
            _userManager = userManager;
            _dbContext = dbContext;
        }

        public IList<Review> Reviews { get; set; }

        public async Task<IActionResult> OnGetAsync()
        {
            var user = await _userManager.GetUserAsync(User);
            if (user == null)
            {
                return NotFound();
            }

            Reviews = await _dbContext.Reviews
                .Include(r => r.Rental)
                .Include(r => r.Reviewer)
                .Where(r => r.Rental.OwnerName.Id == user.Id) 
                .ToListAsync();

            return Page();
        }
    }
}
